file1 = open("readme.dat","w")

Nisomers = int(raw_input("Enter the number of isomers: "))

file1.write("#" + "\n" + "#" + "\n" + str(Nisomers) + "\t" + "|Nisomers" + "\n")

for i in range(1,Nisomers+1):
    file = open("%r.xyz" %i, "r")
#    print(file)
    file_name = ("%r.xyz" %i)
#    print(file_name)

    line = file.read().splitlines()

    line2 = line[1].split()
    energy = line2[0]
#    print(energy)

    file1.write(str(file_name) + "\t" + str(energy) + "\n")

    file.close()

file1.close()
