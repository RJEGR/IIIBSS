def concatenate(file1, file2):
    f1 = open(file1, "r" )
    f2 = open(file2, "r")
    o = open("kareoke.txt", "w")
    content = f1.read() + f2.read()
    o.write(content)
    o.close
    return print("good")
