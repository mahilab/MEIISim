import string

fileheader = ["#pragma once\n",
              "#include <Mahi/Util.hpp>\n",
              "#include <Eigen/Dense>\n",
              "using namespace mahi::util;\n",
              "\n"]

variables = ["\tdouble qe = qs[0];\n",
             "\tdouble qf = qs[1];\n",
             "\tdouble l1 = qs[2];\n",
             "\tdouble l2 = qs[3];\n",
             "\tdouble l3 = qs[4];\n",
             "\tdouble theta1 = qs[5];\n",
             "\tdouble theta2 = qs[6];\n",
             "\tdouble theta3 = qs[7];\n",
             "\tdouble P_p_x = qs[8];\n",
             "\tdouble P_p_y = qs[9];\n",
             "\tdouble P_p_z = qs[10];\n",
             "\tdouble R_p_x = qs[11];\n",
             "\tdouble R_p_y = qs[12];\n",
             "\tdouble R_p_z = qs[13];\n",
             "\tdouble qe_dot = qs[14];\n",
             "\tdouble qf_dot = qs[15];\n",
             "\tdouble l1_dot = qs[16];\n",
             "\tdouble l2_dot = qs[17];\n",
             "\tdouble l3_dot = qs[18];\n",
             "\tdouble theta1_dot = qs[19];\n",
             "\tdouble theta2_dot = qs[20];\n",
             "\tdouble theta3_dot = qs[21];\n",
             "\tdouble V_p_x = qs[22];\n",
             "\tdouble V_p_y = qs[23];\n",
             "\tdouble V_p_z = qs[24];\n",
             "\tdouble w_p_x = qs[25];\n",
             "\tdouble w_p_y = qs[26];\n",
             "\tdouble w_p_z = qs[27];\n",
             "\tdouble sinqe = sin(qe);\n",
             "\tdouble cosqe = cos(qe);\n",
             "\tdouble sinqf = sin(qf);\n",
             "\tdouble cosqf = cos(qf);\n",
             "\tdouble sintheta1 = sin(theta1);\n",
             "\tdouble costheta1 = cos(theta1);\n",
             "\tdouble sintheta2 = sin(theta2);\n",
             "\tdouble costheta2 = cos(theta2);\n",
             "\tdouble sintheta3 = sin(theta3);\n",
             "\tdouble costheta3 = cos(theta3);\n",
             "\tdouble sinR_p_x = sin(R_p_x);\n",
             "\tdouble cosR_p_x = cos(R_p_x);\n",
             "\tdouble sinR_p_y = sin(R_p_y);\n",
             "\tdouble cosR_p_y = cos(R_p_y);\n",
             "\tdouble sinR_p_z = sin(R_p_z);\n",
             "\tdouble cosR_p_z = cos(R_p_z);\n",
             "\n"]

variables_no_dots = ["\tdouble qe = qs[0];\n",
                     "\tdouble qf = qs[1];\n",
                     "\tdouble l1 = qs[2];\n",
                     "\tdouble l2 = qs[3];\n",
                     "\tdouble l3 = qs[4];\n",
                     "\tdouble theta1 = qs[5];\n",
                     "\tdouble theta2 = qs[6];\n",
                     "\tdouble theta3 = qs[7];\n",
                     "\tdouble P_p_x = qs[8];\n",
                     "\tdouble P_p_y = qs[9];\n",
                     "\tdouble P_p_z = qs[10];\n",
                     "\tdouble R_p_x = qs[11];\n",
                     "\tdouble R_p_y = qs[12];\n",
                     "\tdouble R_p_z = qs[13];\n",
                     "\tdouble sinqe = sin(qe);\n",
                     "\tdouble cosqe = cos(qe);\n",
                     "\tdouble sinqf = sin(qf);\n",
                     "\tdouble cosqf = cos(qf);\n",
                     "\tdouble sintheta1 = sin(theta1);\n",
                     "\tdouble costheta1 = cos(theta1);\n",
                     "\tdouble sintheta2 = sin(theta2);\n",
                     "\tdouble costheta2 = cos(theta2);\n",
                     "\tdouble sintheta3 = sin(theta3);\n",
                     "\tdouble costheta3 = cos(theta3);\n",
                     "\tdouble sinR_p_x = sin(R_p_x);\n",
                     "\tdouble cosR_p_x = cos(R_p_x);\n",
                     "\tdouble sinR_p_y = sin(R_p_y);\n",
                     "\tdouble cosR_p_y = cos(R_p_y);\n",
                     "\tdouble sinR_p_z = sin(R_p_z);\n",
                     "\tdouble cosR_p_z = cos(R_p_z);\n",
                     "\n"]

replacements = ["qe", "qf", "theta1", "theta2", "theta3", "R_p_x", "R_p_y", "R_p_z"]

filenames = [    "M",     "V", "G", "psi", "psi_dq", "psi_dq_dt"]

for filename in filenames:
    read_file = open("DynamicEqs/"+filename+".txt")
    content = read_file.read()
    read_file.close()
    
    for replacement in replacements:
        mat_name = filename
        content = string.replace(content,"  T","\t" + mat_name)
        content = string.replace(content,"][",",")
        content = string.replace(content,"[","(")
        content = string.replace(content,"]",")")
        content = string.replace(content,"sin("+replacement+")","sin"+replacement)
        content = string.replace(content,"cos("+replacement+")","cos"+replacement)
        
    write_file = open("../include/" + mat_name + ".hpp","w")
    write_file.writelines(fileheader)
    if mat_name not in ("G", "psi"): 
        write_file.write("inline Eigen::MatrixXd get_" + mat_name + "(const std::vector<double>& qs){\n")
        write_file.write("\tEigen::MatrixXd " + mat_name + " = Eigen::MatrixXd::Zero(14,14); \n\n")

    else: 
        write_file.write("inline Eigen::VectorXd get_" + mat_name + "(const std::vector<double>& qs){\n")
        write_file.write("\tEigen::VectorXd " + mat_name + " = Eigen::VectorXd::Zero(14); \n\n")

    if "V" in mat_name or "dt" in mat_name: write_file.writelines(variables)

    else: write_file.writelines(variables_no_dots)

    write_file.write(content + "\n")
    write_file.write("\treturn " + mat_name + ";\n}")

    write_file.close()