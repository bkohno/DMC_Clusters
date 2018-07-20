file1 = open('timedata_out.dat','r')
file2 = open('timedata_analysis.dat', 'w')

count = 1

for line in file1:
    count+=1
    if count % 2 == 0:
        fields = line.split()
#    print(fields)

        Nw = raw_input("Enter the number of random walkers: ")
#        Cores = raw_input("Enter the number of cores used: ")
        Cores = 32

        field1 = fields[0]
        field2 = fields[1]
        field3 = fields[2]
        field4 = fields[3]

        time = field3[0:8]
#        print(time)
#        print(field3[7])

        if field3[7] == "e":
            time = field3[0:7]

        sections = time.split(":")
        section1 = sections[0]
        section2 = sections[1]
        section3 = sections[2]

        time = float(section1) + float(section2)/60 + float(section3)/3600
#        print(time)

        CPU = field4[0:4]
#        print(CPU)
        if field4[3] == "%":
            CPU = field4[0:3]

        CPU_eff = float(CPU)/(Cores*100)

        file2.write(str(Nw) + "\t" + str(Cores) + "\t" + str(time) + "\t" + str(CPU) + "\t" + str(CPU_eff) + "\n")

#        print(time)

file1.close()
file2.close()
