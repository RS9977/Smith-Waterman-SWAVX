src = SWAVX_TOP.c SWAVX_utils.c SWAVX_256_Func_SeqToSeq_SubMat.c SWAVX_SubMat.c
AVX2_32_D:
	gcc $(src) -mavx2 -lgomp -o SWAVX -D DEBUG
AVX2_32:
	gcc $(src) -mavx2 -lgomp -o SWAVX
clean: 
	rm SWAVX