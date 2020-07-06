# 2020.07.01.1400CDT: Dakai Zhu: Adding a few lines to help test python


t1:
	touch nohup.out
	\rm -rf nohup.out
	nohup python  t1.py all.vcf phenped.txt t1_out.txt &

clean:
	touch t1_out.txt nohup.out
	\rm -rf t1_out.txt nohup.out
