import matplotlib.pyplot as plt
import numpy as np



header2= """
<!DOCTYPE html>
<html>
<head>
<style>
.coff {
    border-spacing: 0px;
    border-collapse:collapse;
    border:1px solid black;
}



.coff table, td, th {text-align: center;}
.coff th, td {border: 2px solid black;}
.coff tr, td {width: 60px;height: 60px}
.coff table {padding: 0 0 0 0;}
.coff th {background-color: #4CAF50;    color: white;}
.coff tr:nth-child(1) {background-color: gray}
.coff td:nth-child(1) {background-color: gray}




.table1 {
    border-spacing: 0px;
    border-collapse:collapse;
    border:1px solid black;
}
.table1 table, td, th {text-align: center;}
.table1 th, td {border: 2px solid black;}
.table1 tr, td {width: 60px;height: 20px}
.table1 table {padding: 0 0 0 0;}
.table1 th {background-color: #4CAF50;    color: white;}
.table1 tr:nth-child(1) {background-color: gray}


.table2 {
    border-spacing: 0px;
    border-collapse:collapse;
    border:1px solid black;
}
.table2 table, td, th {text-align: center;}
.table2 th, td {border: 2px solid black;}
.table2 tr, td {width: 60px;height: 20px}
.table2 table {padding: 0 0 0 0;}
.table2 th {background-color: #4CAF50;    color: white;}
.table2 tr:nth-child(1) {background-color: gray}
.table2 td:nth-child(1) {background-color: gray}


</style>

</head>
<body>
"""


footer2 = """
</body>
"""

def header():
	return header2
def footer():
	return footer2
def cofftohtml(data, header):

	mylist = data.tolist()
	xx = 0
	for x in range(len(mylist)):
		for y in range(len(mylist[0])):
			mylist[x][y] = round(mylist[x][y],2)

	data = mylist[:]
	data.insert(0,header)

	#header.insert(0,"")
	temp = data[:]

	ind = 0
	for xx in temp[0]:
                #print xx
		data[ind].insert(0,header[ind])
		ind += 1


	data[0][0] = ""


        ttt = ""
        for x in data:
                ttt += "<tr>\n"
                for y in x:
                        ttt += "<td>{0}</td>\n".format(y)

                ttt += "</tr>\n"




        table =  "<table class=coff>" + ttt + "</table>"


	return table



def listtohtml(data, tclass = "table1"):

        ttt = ""
        for x in data:
                ttt += "<tr>\n"
                for y in x:
                        ttt += "<td>{0}</td>\n".format(y)
                ttt += "</tr>\n"
        table =  "<table class=%s>"%tclass + ttt + "</table>"
	return table


#header= ["a1","a2","a3","a4","a5","a6","a7","a8","a9",]
#tohtml(data,header)
















