src = SWAVX_utils.c SWAVX_SubMat.c
AVX2_32_D:
	gcc -O3 $(src) SWAVX_TOP.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -lgomp -o SWAVX -D DEBUG
AVX2_32:
	gcc -O3 $(src) SWAVX_TOP.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -lgomp -o SWAVX

AVX2_32_D_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -lgomp -o SWAVX -D DEBUG
AVX2_32_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -lgomp -o SWAVX

clean: 
	rm SWAVX