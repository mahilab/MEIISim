import string

fileheader = ["#pragma once\n","#include <Mahi/Util.hpp>\n","using namespace mahi::util;\n","\n"]

variables = ["\tdouble l1 = qs[0];\n","\tdouble l2 = qs[1];\n","\tdouble l3 = qs[2];\n","\tdouble theta1 = qs[3];\n","\tdouble theta2 = qs[4];\n","\tdouble theta3 = qs[5];\n","\tdouble l1_dot = qs[6];\n","\tdouble l2_dot = qs[7];\n","\tdouble l3_dot = qs[8];\n","\tdouble theta1_dot = qs[9];\n","\tdouble theta2_dot = qs[10];\n","\tdouble theta3_dot = qs[11];\n","\tdouble sintheta1 = sin(theta1);\n","\tdouble costheta1 = cos(theta1);\n","\tdouble sintheta2 = sin(theta2);\n","\tdouble costheta2 = cos(theta2);\n","\tdouble sintheta3 = sin(theta3);\n","\tdouble costheta3 = cos(theta3);\n","\n"]

variables_no_dots = ["\tdouble l1 = qs[0];\n","\tdouble l2 = qs[1];\n","\tdouble l3 = qs[2];\n","\tdouble theta1 = qs[3];\n","\tdouble theta2 = qs[4];\n","\tdouble theta3 = qs[5];\n","\tdouble sintheta1 = sin(theta1);\n","\tdouble costheta1 = cos(theta1);\n","\tdouble sintheta2 = sin(theta2);\n","\tdouble costheta2 = cos(theta2);\n","\tdouble sintheta3 = sin(theta3);\n","\tdouble costheta3 = cos(theta3);\n","\n"]

replacements = ["theta1", "theta2", "theta3"]


filenames = ["Mpar11", "Mpar12", "Mpar13", "Mpar21", "Mpar22", "Mpar23", "Mpar31", "Mpar32", "Mpar33", "Vpar11", "Vpar12", "Vpar13", "Vpar21", "Vpar22", "Vpar23", "Vpar31", "Vpar32", "Vpar33", "Gpar1", "Gpar2", "Gpar3"]

for filename in filenames:
    read_file = open("DynamicEqs/"+filename+".txt")
    content = read_file.read()
    read_file.close()
    
    for replacement in replacements:
        mat_name = string.replace(filename,"par","")
        content = string.replace(content,"  t0","\tdouble " + mat_name)
        content = string.replace(content,"sin("+replacement+")","sin"+replacement)
        content = string.replace(content,"cos("+replacement+")","cos"+replacement)
        
    write_file = open("../include/" + mat_name + ".hpp","w")
    write_file.writelines(fileheader)
    write_file.write("double get_" + mat_name + "(std::vector<double> qs){\n")
    if "V" in mat_name: write_file.writelines(variables)
    else: write_file.writelines(variables_no_dots)
    write_file.write(content + "\n")
    write_file.write("\treturn " + mat_name + ";\n}")