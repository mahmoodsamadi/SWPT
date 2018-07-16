import arcpy,os,sys,math
from arcpy import env
from arcpy.sa import *
arcpy.env.overwriteOutput = True
	
	
inRasters = arcpy.GetParameterAsText(0)
workspace = arcpy.GetParameterAsText(1)
convalue = 1500



#############################################################
def doit(workspace,inRaster):

	dec= arcpy.Describe(inRaster)
	#if dec.spatialReference.linearUnitName != 'Meter':
	#	arcpy.AddError("raster({0}) linear Unit is not Meter;=>{1} .".format(dec.baseName, str(dec.spatialReference.linearUnitName)))
	#	sys.exit(0)

	workspace = workspace + "/" + dec.baseName
	if not os.path.exists(workspace):
		os.makedirs(workspace)

	dscRaster = arcpy.Describe(inRaster)
	for child in dscRaster.children:
		cellSize = child.meanCellHeight
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic = {}
	vardic["name"] = dec.baseName
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	item = Con( Raster(inRaster) > 0, 1, 1)
	item.save(workspace + "/item")
	#item = Divide(temprr, temprr)
	
	out_basin_tmp = workspace + "/basin_tmp.shp"
	out_basin = workspace + "/basin.shp"
	
	dir_path = workspace + "\\f_DIRECTION"
	acc_path = workspace + "\\f_Accumul"

	
	###################################
	def addfield(inFeatures, fieldname):
		arcpy.AddField_management(inFeatures, fieldname, "DOUBLE")
	
	def sss(inFeatures, fieldName1,what):
		addfield(inFeatures, fieldName1)
		arcpy.CalculateField_management(inFeatures, fieldName1,what,"PYTHON_9.3")
	
	def field2list(infc,field):
		list = []
		for row in arcpy.da.SearchCursor(infc, [field]):
			list.append(row[0])
		return list
	
	
	
	#####################################
	
	demfill = Fill(inRaster)

	outFlowDirection = FlowDirection(demfill, "FORCE", "")
	outFlowDirection.save(dir_path)
	
	outFlowAccumulation = FlowAccumulation(outFlowDirection, "", "INTEGER")
	outFlowAccumulation.save(acc_path)
	
	
	
	######## get TWI,SPI,sti

	mslope = arcpy.sa.Slope(demfill)
	mslope = Con( mslope > 1, mslope, 1 )

	slope_radians = mslope * math.pi/180.0
	slope_radians = Con( slope_radians > 0, slope_radians, 0.001 )

	sca = ((outFlowAccumulation + 1) * cellSize * cellSize)
	TWI = arcpy.sa.Ln(sca / (arcpy.sa.Tan(slope_radians)))
	TWI = Con( TWI > 0, TWI, 0.01 )
	TWI.save(workspace + "/twi")

	twi_mean = arcpy.GetRasterProperties_management(TWI,"MEAN").getOutput(0)
	twi_mean = float(twi_mean.replace(",","."))
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["TWI"] = twi_mean
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	#Stream power index (SI)*
	SPI = arcpy.sa.Ln(sca * (arcpy.sa.Tan(slope_radians)))
	SPI = Con( SPI > 0, SPI, 0.01 )
	SPI.save(workspace + "/spi")

	spi_mean = arcpy.GetRasterProperties_management(SPI ,"MEAN").getOutput(0)
	spi_mean = float(spi_mean.replace(",","."))
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["SPI"] = spi_mean 
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


	#Sediment transport index
	STI = (0.4 + 1.0) * ((sca/22.13)**0.4) * (arcpy.sa.Sin(slope_radians/0.0896)**0.0896)
	STI.save(workspace + "/sti")

	sti_mean = arcpy.GetRasterProperties_management(STI ,"MEAN").getOutput(0)
	sti_mean = float(sti_mean.replace(",","."))
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["STI"] = sti_mean 
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



	# getting outlet x,y
	#######################################
	
	acc_max = arcpy.GetRasterProperties_management(acc_path,"MAXIMUM").getOutput(0)
	acc_max = acc_max.replace(",",".")
	
	accmaxpoint = workspace + "/accmax_point.shp"
	outCon = Con(Raster(acc_path) == float(acc_max),1, 0)
	outSetNull = SetNull(outCon, 1, "VALUE <> 1")
	field = "VALUE"
	arcpy.RasterToPoint_conversion(outSetNull, accmaxpoint, field)
	
	sss(accmaxpoint, "x","!shape.centroid.X!")
	sss(accmaxpoint, "y","!shape.centroid.Y!")
	
	x = field2list(accmaxpoint,"x")
	y = field2list(accmaxpoint,"y")
	
	out_x,out_y = float(x[0]),float(y[0])
	
	########################################
	arcpy.RasterToPolygon_conversion(item, out_basin_tmp, "NO_SIMPLIFY","VALUE")
	arcpy.Dissolve_management(out_basin_tmp, out_basin, "", "", "MULTI_PART", "DISSOLVE_LINES")

	

	mlist =[]
	def polygons(inFC):
		for row in arcpy.da.SearchCursor(inFC, ["OID@", "SHAPE@"]):
			for part in row[1]:
				for vertex in part:
					mlist.append([float(vertex.X), float(vertex.Y)])
	arcpy.AddMessage(out_basin)
	polygons(out_basin)
	#arcpy.AddMessage(mlist)

	
	minVal = arcpy.GetRasterProperties_management(inRaster,"MINIMUM").getOutput(0)
	minVal = minVal.replace(",",".")
	maxVal = arcpy.GetRasterProperties_management(inRaster,"MAXIMUM").getOutput(0)
	maxVal = maxVal.replace(",",".")
	
	
	#minVal
	outminpoint2 = workspace + "/min_h.shp"
	outminpoint = workspace + "/min_h_single.shp"

	outCon = Con(Raster(inRaster) == float(minVal),1, 0)
	outSetNull = SetNull(outCon, 1, "VALUE <> 1")


	field = "VALUE"
	arcpy.RasterToPoint_conversion(outSetNull, outminpoint2, field)
	arcpy.Select_analysis(outminpoint2, outminpoint,'"FID" = 0')


	sss(outminpoint, "x","!shape.centroid.X!")
	sss(outminpoint, "x","!shape.centroid.Y!")
	
	
	#maxVal
	outmaxpoint2 = workspace + "/max_h.shp"
	outmaxpoint = workspace + "/max_h_single.shp"

	outCon = Con(Raster(inRaster) == float(maxVal),1, 0)
	outSetNull = SetNull(outCon, 1, "VALUE <> 1")
	field = "VALUE"
	arcpy.RasterToPoint_conversion(outSetNull, outmaxpoint2, field)
	arcpy.Select_analysis(outmaxpoint2, outmaxpoint,'"FID" = 0')
	
	
	# River ##################################

	minacc = arcpy.GetRasterProperties_management(outFlowAccumulation,"MINIMUM").getOutput(0)
	minacc = int(minacc.replace(",","."))
	maxacc = arcpy.GetRasterProperties_management(outFlowAccumulation,"MAXIMUM").getOutput(0)
	maxacc = int(maxacc.replace(",","."))






	#accminVal
	outminpoint2 = workspace + "/min_acc.shp"
	outminpoint = workspace + "/min_acc_single.shp"

	outCon = Con(outFlowAccumulation == float(minacc),1, 0)
	outSetNull = SetNull(outCon, 1, "VALUE <> 1")


	field = "VALUE"
	arcpy.RasterToPoint_conversion(outSetNull, outminpoint2, field)
	arcpy.Select_analysis(outminpoint2, outminpoint,'"FID" = 0')


	sss(outminpoint, "x","!shape.centroid.X!")
	sss(outminpoint, "x","!shape.centroid.Y!")
	
	
	#accmaxVal
	outmaxpoint2 = workspace + "/max_acc.shp"
	outmaxpoint = workspace + "/max_acc_single.shp"

	outCon = Con(outFlowAccumulation == float(maxacc),1, 0)
	outSetNull = SetNull(outCon, 1, "VALUE <> 1")
	field = "VALUE"
	arcpy.RasterToPoint_conversion(outSetNull, outmaxpoint2, field)
	arcpy.Select_analysis(outmaxpoint2, outmaxpoint,'"FID" = 0')
	








	percentage = 0.5
	convalue = (percentage * (maxacc - minacc) / 100) + minacc
	arcpy.AddMessage("@"*33)
	

	arcpy.AddMessage((100 * (maxacc - minacc) / 100) + minacc)


	vall = outFlowAccumulation > int(convalue)
	streams = Con(vall , 1, "")
	streams.save(workspace + "\\river")
	
	outStreamOrder = StreamOrder(streams, outFlowDirection, "STRAHLER")
	outStreamOrder.save(workspace + "\\StrOrdSTRAH")
	strorder = workspace + "\\StrOrdSTRAH.shp"
	StreamToFeature(outStreamOrder, outFlowDirection, strorder, "NO_SIMPLIFY")

	################################################################################################################################################################################
	#for row in arcpy.da.SearchCursor(inFC, ["TO_NODE", "SHAPEeeeee@","GRID_CODE"]):
	#	print 1111
	###
	def line(inFC,workspace,out_spliterpoints):
		mlist = []
		for row in arcpy.da.SearchCursor(inFC, ["TO_NODE", "SHAPE@","GRID_CODE"]):
			mm = row[0]
			order = row[2]
			#xy_start = [row[1].positionAlongLine(0.0,True).firstPoint.X, row[1].positionAlongLine(0.0,True).firstPoint.Y]
			#mlist.append([mm,xy_start])
			xy_end = [row[1].positionAlongLine(1.0,True).firstPoint.X, row[1].positionAlongLine(1.0,True).firstPoint.Y]
			mlist.append([str(mm)+str(order),xy_end,order])
		nods = mlist
	
		mmlist = []
		points = []
		for nod in nods:
			if nod[0] in mmlist:
				points.append(nod[1])
			mmlist.append(nod[0])
	
		outPath, outName = os.path.split(out_spliterpoints)
	
		arcpy.CreateFeatureclass_management(outPath, outName, "Multipoint", "", "DISABLED", "DISABLED", strorder)
		print len(points)
		cursor = arcpy.da.InsertCursor(out_spliterpoints, ("SHAPE@"))
		for xy in points:
			centroidx = xy[0]
			centroidy = xy[1]
			print centroidx,centroidy
			array = arcpy.Array(arcpy.Point(centroidx,centroidy))
			Multipoint = arcpy.Multipoint(array)
			cursor.insertRow((Multipoint,))
		del cursor
	
	##############
	
	strordstrah_Dissolve2 = workspace + "\\str_diss.shp"
	arcpy.Dissolve_management(strorder, strordstrah_Dissolve2, "GRID_CODE", "", "MULTI_PART", "DISSOLVE_LINES")
	
	all_orders = field2list(strordstrah_Dissolve2, "GRID_CODE")
	for orders in all_orders:
		arcpy.AddMessage("#"*33)
		arcpy.AddMessage(orders)
	
		ppath = workspace + "\\ord_%s.shp"%orders
		arcpy.Select_analysis(strorder, ppath, "\"GRID_CODE\" =%s"%orders)
	
		out_spliterpoints = workspace + "\\spliterpoint%s.shp"%orders
	
		line(ppath,workspace,out_spliterpoints)
	
		order_diss = workspace + "\\ord_%s_diss.shp"%orders
		arcpy.Dissolve_management(ppath, order_diss, "GRID_CODE", "", "MULTI_PART", "DISSOLVE_LINES")
	
		mystrahler = workspace + "\\ord_spl_%s.shp"%orders
		arcpy.SplitLineAtPoint_management(order_diss, out_spliterpoints, mystrahler, "20 Meters")
	
	
	
	order_info = {}
	
	for orders in all_orders:
		mystrahler = workspace + "\\ord_spl_%s.shp"%orders
		sss(mystrahler, "length","!shape.length@METERS!")
		x = field2list(mystrahler, "length")
		length = 0
		for segment in x:
			length += float(segment)
		
		count = len(x)
		order_info[orders] = [count,length]
	
	arcpy.AddMessage(order_info)
	
	

	
	####################################
	####################################
	# Nu ########################
	Nu = 0
	for x in order_info:
		Nu += order_info[x][0]
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["Nu"] = Nu
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	
	# Lu ########################
	Lu = 0
	for x in order_info:
		Lu += order_info[x][1]
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["Lu"] = Lu
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	
	
	# Rb ################################
	rb= 0
	for x in range(1,len(order_info)):
		rb += float(order_info[x][0])/float(order_info[x+1][0])
	Rb =  rb / len(order_info) - 1
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["Rb"] = Rb
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	
	# Rl ################################
	rl= 0
	for x in range(1,len(order_info)):
		rl += float(order_info[x][1])/float(order_info[x+1][1])
	Rl =  rl / len(order_info) - 1
	
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["Rl"] = Rl
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	
	# Lb ###############################
	
	mydistance = 0
	for xy in mlist:
		dis = ((xy[0] - out_x)**2 + (xy[1] - out_y)**2) ** 0.5
		if dis > mydistance:
			mydistance = dis
			far_x,far_y = xy[0],xy[1]
	
	
	outFC = workspace + "/" + "ffffffff.shp"
	outPath, outName = os.path.split(outFC)
	arcpy.CreateFeatureclass_management(outPath, outName, "Polyline", "", "DISABLED", "DISABLED", out_basin)
	array = arcpy.Array()
	
	array.add(arcpy.Point(far_x,far_y))
	array.add(arcpy.Point(out_x,out_y))
	
	Polyline = arcpy.Polyline(array)
	cursor = arcpy.da.InsertCursor(outFC, ("SHAPE@"))
	cursor.insertRow((Polyline,))
	del cursor
	
	sss(outFC, "length","!shape.length@METERS!")
	x = field2list(outFC, "length")
	Lb = 0
	for segment in x:
		Lb += float(segment)
	
	arcpy.AddMessage("#Lb"*33)
	arcpy.AddMessage(Lb)
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["Lb"] = Lb
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	
	# Lm ##############################################
	field = "VALUE"
	outFlowLength = FlowLength(outFlowDirection, "DOWNSTREAM", "")
	mmmVal = arcpy.GetRasterProperties_management(outFlowLength,"MAXIMUM").getOutput(0)
	mmmVal = mmmVal.replace(",",".")
	outCon = Con(outFlowLength == float(mmmVal),1, 0)
	outSetNull = SetNull(outCon, 1, "VALUE <> 1")
	
	outPoint = workspace + "/" + "l_f_p_point.shp"
	outstream = workspace + "/" + "l_f_p.shp"
	
	arcpy.RasterToPoint_conversion(outSetNull, outPoint, field)
	outCostPath = CostPath(outPoint, outFlowDirection, outFlowDirection, "EACH_CELL")
	StreamToFeature(outCostPath, outFlowDirection,outstream, "NO_SIMPLIFY")
	
	
	
	sss(outstream, "length","!shape.length@METERS!")
	x = field2list(outstream, "length")
	Lm = 0
	for segment in x:
		Lm += float(segment)
	arcpy.AddMessage("#Lm"*33)
	arcpy.AddMessage(Lm)
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["Lm"] = Lm
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	
	# A #################################################
	sss(out_basin, "area","!shape.area@SQUAREMETERS!")
	x = field2list(out_basin, "area")
	A = 0
	for segment in x:
		A += float(segment)
	arcpy.AddMessage("#A"*33)
	arcpy.AddMessage(A)
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["A"] = A
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	# P #################################################
	sss(out_basin, "perimeter","!shape.length@METERS!")
	x = field2list(out_basin, "perimeter")
	P = 0
	for segment in x:
		P += float(segment)
	arcpy.AddMessage("#P"*33)
	arcpy.AddMessage(P)
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["P"] = P
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	
	# H #################################################
	H = float(maxVal) - float(minVal)
	arcpy.AddMessage("#H"*33)
	arcpy.AddMessage(H)
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["H"] = H
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	# Cc #######################
	sss(out_basin, "Cc","0.2821 * (float(!perimeter!) / (float(!area!) ** 0.5))")
	x = field2list(out_basin, "Cc")
	Cc = 0
	for segment in x:
		Cc += float(segment)
	arcpy.AddMessage("#Cc"*33)
	arcpy.AddMessage(Cc)
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["Cc"] = Cc
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	# Rc ########################
	sss(out_basin, "Rc","12.56 * (float(!area!) / (float(!perimeter!) ** 2))")
	x = field2list(out_basin, "Rc")
	Rc = 0
	for segment in x:
		Rc += float(segment)
	arcpy.AddMessage("#Rc"*33)
	arcpy.AddMessage(Rc)
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["Rc"] =Rc
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	# Rf #########################
	sss(out_basin, "Rf","(float(!area!) / (float(%s) ** 2))" % Lb)
	x = field2list(out_basin, "Rf")
	Rf = 0
	for segment in x:
		Rf += float(segment)
	arcpy.AddMessage("#Rf"*33)
	arcpy.AddMessage(Rf)
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["Rf"] =Rf
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	# Re ##########################
	sss(out_basin, "Re","2 * (((float(!area!)/3.1415) ** 0.5)/ (float(%s)))" % Lb)
	x = field2list(out_basin, "Re")
	Re = 0
	for segment in x:
		Re += float(segment)
	arcpy.AddMessage("#Re"*33)
	arcpy.AddMessage(Re)
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["Re"] =Re
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	
	# Fs ###############################
	Fs = Nu / A
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["Fs"] =Fs
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	# D ###############################
	D = Lu / A
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["D"] =D
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	# Lg ###############################
	Lg = 1 / (D * 2)
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["Lg"] =Lg
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	# Lsm ##############################
	Lsm = Lu / Nu
	
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["Lsm"] =Lsm
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	# Rt ##############################
	Rt = Nu/P
	
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["Rt"] =Rt
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	# Rh ##############################
	Rh = H/Lb
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["Rh"] =Rh
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	# Rn ##############################
	Rn = H * D
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["Rn"] =Rn
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	# C ##############################
	C = 1/D
	
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	#@
	vardic["C"] =C
	#@
	#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	
	arcpy.AddMessage("@"*33)
	for vvv in vardic:
		arcpy.AddMessage("%s,%s"%(vvv,vardic[vvv]))
	
	return vardic
	
##################################################################################################################################################




data =[]

inRasters = inRasters.split(";")
for inRaster in inRasters:
	data.append(doit(workspace,inRaster))
arcpy.AddMessage("data is:")
arcpy.AddMessage(data)
arcpy.AddMessage("#"* 33)

all_vars =["name","Rb","Rf","Re","Rc","D","Rt","Cc","C","Lg","Lb","Lm","Lu","Rl","Rn","Rh","Nu","Fs","A","H","P","Lsm","TWI", "SPI", "STI"]
#this is to sort


myindices = ["Fs","Rb","Rf","Re","Rc","D","Rt","Cc","C","TWI", "SPI", "STI"]
mat =[]
for xx in myindices:
	temp = []
	for dem in data:
		myvar = dem[xx]
		#if dem[xx] < 0.00005: myvar = 0
		temp.append(myvar)
		
	mat.append(temp)

arcpy.AddMessage("mat is:")
arcpy.AddMessage(mat)
arcpy.AddMessage("#"* 33)

import numpy as np
data_cof = np.corrcoef(mat)

arcpy.AddMessage(data_cof)

arcpy.AddMessage("#"* 33)





###########################################
sum_of_var = list(np.sum(data_cof, axis=0))
sum_of_all = float(np.sum(data_cof))
sum_dic = dict(zip(myindices, sum_of_var))

arcpy.AddMessage("#"* 33)
arcpy.AddMessage("sumall")
arcpy.AddMessage(sum_dic)
arcpy.AddMessage("#"* 33)

result = []
for dem in data:
		temp_dic = {}
		Prioritization = 0.0
		temp_dic["name"] = dem["name"]
		for var in myindices:
			x = sum_dic[var] / sum_of_all * dem[var]
			arcpy.AddMessage(sum_dic[var])
			arcpy.AddMessage(sum_of_all)
			arcpy.AddMessage(dem[var])
			Prioritization += x
		temp_dic["Prioritization"] = Prioritization
		result.append(temp_dic)

arcpy.AddMessage("#"* 33)
arcpy.AddMessage("#"* 33)
arcpy.AddMessage("#"* 33)
arcpy.AddMessage(result)
arcpy.AddMessage("#"* 33)
arcpy.AddMessage("#"* 33)
arcpy.AddMessage("#"* 33)
################
arcpy.AddMessage(data)
f = open(workspace + "/report.txt","w")
f.write(",".join(all_vars))
f.write("\n")
for dem in data:
	for ss in all_vars:
		f.write(str(dem[ss]))
		f.write(",")
	f.write("\n")
f.close()















###########################
#coff table


import matplotlib.pyplot as plt


def show_values(pc, fmt="%.2f", **kw):
    from itertools import izip
    pc.update_scalarmappable()
    ax = pc.get_axes()
    for p, color, value in izip(pc.get_paths(), pc.get_facecolors(), pc.get_array()):
        x, y = p.vertices[:-2, :].mean(0)
        if np.all(color[:3] > 0.5):
            color = (0.0, 0.0, 0.0)
        else:
            color = (1.0, 1.0, 1.0)
        ax.text(x, y, fmt % value, ha="center", va="center", color=color, **kw)


def cm2inch(*tupl):
    inch = 2.54
    if type(tupl[0]) == tuple:
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


def heatmap(AUC, title, xlabel, ylabel, xticklabels, yticklabels, figure_width=40, figure_height=20, correct_orientation=False, cmap='RdBu'):

    # Plot it out
    fig, ax = plt.subplots()    
    #c = ax.pcolor(AUC, edgecolors='k', linestyle= 'dashed', linewidths=0.2, cmap='RdBu', vmin=0.0, vmax=1.0)
    c = ax.pcolor(AUC, edgecolors='k', linestyle= 'dashed', linewidths=0.2, cmap=cmap)

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(AUC.shape[0]) + 0.5, minor=False)
    ax.set_xticks(np.arange(AUC.shape[1]) + 0.5, minor=False)

    # set tick labels
    #ax.set_xticklabels(np.arange(1,AUC.shape[1]+1), minor=False)
    ax.set_xticklabels(xticklabels, minor=False)
    ax.set_yticklabels(yticklabels, minor=False)




    # set title and x/y labels
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)      

    # Remove last blank column
    plt.xlim( (0, AUC.shape[1]) )

    # Turn off all the ticks
    ax = plt.gca()    
    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    # Add color bar
    plt.colorbar(c)

    # Add text in each cell 
    show_values(c)

    # Proper orientation (origin at the top left instead of bottom left)
    if correct_orientation:
        ax.invert_yaxis()
        ax.xaxis.tick_top()       

    # resize 
    fig = plt.gcf()
    #fig.set_size_inches(cm2inch(40, 20))
    #fig.set_size_inches(cm2inch(40*4, 20*4))
    fig.set_size_inches(cm2inch(figure_width, figure_height))



def plot_classification_report(data,mylabel=[], title='', cmap='RdBu'):
    xlabel = ''
    ylabel = ''
    xticklabels = mylabel
    yticklabels = mylabel
    figure_width = 25
    figure_height = len(yticklabels) + 7
    correct_orientation = False
    correct_orientation = True
    heatmap(data, title, xlabel, ylabel, xticklabels, yticklabels, figure_width, figure_height, correct_orientation, cmap=cmap)


plot_classification_report(data_cof,mylabel=myindices)
plt.savefig(workspace + "/report.png",mylabel=myindices, dpi=200, format='png', bbox_inches='tight')
plt.close()



##########################################################

import tohtml
html = tohtml.header()

myindices222 = ["name","Fs","Rb","Rf","Re","Rc","D","Rt","Cc","C","TWI", "SPI", "STI"]
mat =[]
for xx in myindices222:
	temp = [xx]
	for dem in data:
		myvar = dem[xx]
		temp.append(myvar)
	mat.append(temp)
html2 = tohtml.listtohtml(mat, tclass = "table1")
html = html + html2 + "<br>"



html2 = tohtml.cofftohtml(data_cof,myindices)
html = html + html2 + "<br>"



sss =[["Name",'Prioritization']]
for x in result:
	temp = [x['name'],x['Prioritization']]
	sss.append(temp)
html2 = tohtml.listtohtml(sss, tclass = "table1")
html = html + html2 + "<br>"


f = open(workspace + "/report.html","w")
f.write(html)
f.close()


