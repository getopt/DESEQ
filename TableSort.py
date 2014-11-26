import sys

FileName=sys.argv[1]
Quary=sys.argv[2]
OutputContain=sys.argv[3]
OutputNotContain=sys.argv[4]

original=open(FileName,'r')
output_with_str=open(OutputContain,'w')
output_wo_str=open(OutputNotContain,'w')

for line in original:
    if line.find(Quary)==-1:
        output_with_str.write(line)
    else:
        output_wo_str.write(line)
print("done")
