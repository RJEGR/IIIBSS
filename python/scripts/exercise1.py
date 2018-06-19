def UPERCONVERT(filename):
    i = open(filename, "r")
    content = i.read()
    o = open("SUPERtotalEclipse.txt", "w") # create out file
    o.write( content.replace("heart", "moon").upper() )   # write the file by pipe the replace and upper function through your data content
    o.close() # close it :)
    i.close()     
    return print("bye")
