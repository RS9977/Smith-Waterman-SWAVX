src = SWAVX_utils.c SWAVX_SubMat.c
AVX2_32_D_SM:
	gcc -O3 $(src) SWAVX_TOP.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -o SWAVX -D DEBUG -D SUBMAT
AVX2_32_SM:
	gcc -O3 $(src) SWAVX_TOP.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -o SWAVX -D SUBMAT

AVX2_32_D_SM_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -o SWAVX -D DEBUG -D SUBMAT
AVX2_32_SM_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -o SWAVX -D SUBMAT

AVX2_32_D:
	gcc -O3 $(src) SWAVX_TOP.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -o SWAVX -D DEBUG
AVX2_32:
	gcc -O3 $(src) SWAVX_TOP.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -o SWAVX

AVX2_32_D_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -o SWAVX -D DEBUG
AVX2_32_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -o SWAVX


AVX512_32_D_SM:
	gcc -O3 $(src) SWAVX_TOP.c SWAVX_512_Func_SeqToSeq_SubMat.c -mavx512f -o SWAVX -D DEBUG -D SUBMAT -D B512
AVX512_32_SM:
	gcc -O3 $(src) SWAVX_TOP.c SWAVX_512_Func_SeqToSeq_SubMat.c -mavx512f -o SWAVX -D SUBMAT -D B512

AVX512_32_D_SM_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets.c SWAVX_512_Func_SeqToSeq_SubMat.c -mavx512f -o SWAVX -D DEBUG -D SUBMAT -D B512
AVX512_32_SM_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets.c SWAVX_512_Func_SeqToSeq_SubMat.c -mavx512f -o SWAVX -D SUBMAT -D B512

AVX512_32_D:
	gcc -O3 $(src) SWAVX_TOP.c SWAVX_512_Func_SeqToSeq_SubMat.c -mavx512f -o SWAVX -D DEBUG -D B512
AVX512_32:
	gcc -O3 $(src) SWAVX_TOP.c SWAVX_512_Func_SeqToSeq_SubMat.c -mavx512f -o SWAVX -D B512

AVX512_32_D_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets.c SWAVX_512_Func_SeqToSeq_SubMat.c -mavx512f -o SWAVX -D DEBUG -D B512
AVX512_32_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets.c SWAVX_512_Func_SeqToSeq_SubMat.c -mavx512f -o SWAVX -D B512





AVX512_32_D_SM_MT_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets_MultiThread.c SWAVX_512_Func_SeqToSeq_SubMat.c -mavx512f -lpthread -o SWAVX -D DEBUG -D SUBMAT -D B512
AVX512_32_SM_MT_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets_MultiThread.c SWAVX_512_Func_SeqToSeq_SubMat.c -mavx512f -lpthread -o SWAVX -D SUBMAT -D B512

AVX512_32_D_MT_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets_MultiThread.c SWAVX_512_Func_SeqToSeq_SubMat.c -mavx512f -lpthread -o SWAVX -D DEBUG -D B512
AVX512_32_MT_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets_MultiThread.c SWAVX_512_Func_SeqToSeq_SubMat.c -mavx512f -lpthread -o SWAVX -D B512

AVX2_32_D_SM_MT_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets_MultiThread.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -lpthread -o SWAVX -D DEBUG -D SUBMAT
AVX2_32_SM_MT_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets_MultiThread.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -lpthread -o SWAVX -D SUBMAT

AVX2_32_D_MT_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets_MultiThread.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -lpthread -o SWAVX -D DEBUG
AVX2_32_MT_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets_MultiThread.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -lpthread -o SWAVX

clean: 
	rm SWAVX