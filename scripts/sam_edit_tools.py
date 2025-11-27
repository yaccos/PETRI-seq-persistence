import os

def uniform_sam(sample):

	input_sam = open(sample)
	#input_sam = open(sample + '.sam')
	output = open(sample+'_uniform_XT.sam','w')

	for line in input_sam:
		if("XT:A:U" in line):
			output.write(line.replace("XT:A:U", "XT:A:N"))
		elif("XT:A:R" in line):
			output.write(line.replace("XT:A:R", "XT:A:N"))
		else:
			output.write(line)
	output.close()


def no_xt(sample):

	input_sam = open(sample + '.sam')
	output = open(sample+'_no_xt.sam','w')
	for line in input_sam:
		if("XT:" in line):
			output.write(line.replace("XT:", "XN:"))
		else:
			output.write(line)
	output.close()

def no_xt_new(sample_in,sample_out):
	input_sam = open(sample_in)
	output = open(sample_out,'w')
	for line in input_sam:
		if("XT:" in line):
			output.write(line.replace("XT:", "XN:"))
		else:
			output.write(line)
	output.close()
