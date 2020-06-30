import string

fileheader = ["#pragma once\n",
              "#include <Mahi/Util.hpp>\n",
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
             "\tdouble qe_dot = qs[8];\n",
             "\tdouble qf_dot = qs[9];\n",
             "\tdouble l1_dot = qs[10];\n",
             "\tdouble l2_dot = qs[11];\n",
             "\tdouble l3_dot = qs[12];\n",
             "\tdouble theta1_dot = qs[13];\n",
             "\tdouble theta2_dot = qs[14];\n",
             "\tdouble theta3_dot = qs[15];\n",
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
             "\n"]

variables_no_dots = ["\tdouble qe = qs[0];\n",
                     "\tdouble qf = qs[1];\n",
                     "\tdouble l1 = qs[2];\n",
                     "\tdouble l2 = qs[3];\n",
                     "\tdouble l3 = qs[4];\n",
                     "\tdouble theta1 = qs[5];\n",
                     "\tdouble theta2 = qs[6];\n",
                     "\tdouble theta3 = qs[7];\n",
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
                     "\n"]

replacements = ["qe", "qf", "theta1", "theta2", "theta3"]

filenames = ["M11", "M12", "M13", "M14", "M15", "M16", "M17", "M18",
             "M21", "M22", "M23", "M24", "M25", "M26", "M27", "M28",
             "M31", "M32", "M33", "M34", "M35", "M36", "M37", "M38",
             "M41", "M42", "M43", "M44", "M45", "M46", "M47", "M48",
             "M51", "M52", "M53", "M54", "M55", "M56", "M57", "M58",
             "M61", "M62", "M63", "M64", "M65", "M66", "M67", "M68",
             "M71", "M72", "M73", "M74", "M75", "M76", "M77", "M78",
             "M81", "M82", "M83", "M84", "M85", "M86", "M87", "M88",
             "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18",
             "V21", "V22", "V23", "V24", "V25", "V26", "V27", "V28",
             "V31", "V32", "V33", "V34", "V35", "V36", "V37", "V38",
             "V41", "V42", "V43", "V44", "V45", "V46", "V47", "V48",
             "V51", "V52", "V53", "V54", "V55", "V56", "V57", "V58",
             "V61", "V62", "V63", "V64", "V65", "V66", "V67", "V68",
             "V71", "V72", "V73", "V74", "V75", "V76", "V77", "V78",
             "V81", "V82", "V83", "V84", "V85", "V86", "V87", "V88",
             "G1",   "G2",  "G3",  "G4",  "G5",  "G6",  "G7",  "G8",
             "rho", "rhodt"]

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
    write_file.write("inline double get_" + mat_name + "(const std::vector<double>& qs){\n")
    if "V" in mat_name: write_file.writelines(variables)
    else: write_file.writelines(variables_no_dots)
    write_file.write(content + "\n")
    write_file.write("\treturn " + mat_name + ";\n}")