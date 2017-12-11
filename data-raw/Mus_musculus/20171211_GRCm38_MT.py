"""
Created on Sat Oct 14 17:01:47 2017

@author: Zhe
"""

def get_MT(file, out):
    with open(file, mode = "r") as f:
        with open(out, mode = "w") as o:
	    flag = 0
            for i in f:
                if i.startswith(">MT"):
		    flag = 1
                    o.write(i)
		elif flag == 1:
		    o.write(i)
		

def main():
    file = "/restricted/projectnb/cbmhive/references/Mouse/GRCm38/raw/GCA_000001635.8_GRCm38.p6_genomic_modified.fna"
    out = "./GCA_000001635.8_GRCm38.p6_genomic_modified_MT.fna"
    get_MT(file, out)

if __name__ == "__main__":
    main()
