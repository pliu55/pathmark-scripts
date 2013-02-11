import os, re
import mData

assert(os.path.exists("first.props"))
assert(os.path.exists("include.features"))

features = mData.rList("include.features")

f = open("first.props", "r")
o = open("all.props", "w")
for line in f:
    if line.startswith("nodeCustomGraphics1.default-Node\ Custom\ Graphics\ 1-Discrete\ Mapper.mapping.map"):
        pline = re.split(",", line.rstrip("\r\n"))
        print pline
        imgCounter = int(pline[1])
        for feature in features:
            reline = re.sub(features[0], feature, pline[0]) + "," + str(imgCounter) + "," + re.sub(re.sub("/", "_", features[0]), re.sub("/", "_", feature), pline[2]) + "," + pline[3]
            imgCounter += 1
            o.write(reline+"\n")
    else:
        o.write(line)
f.close()
o.close()
