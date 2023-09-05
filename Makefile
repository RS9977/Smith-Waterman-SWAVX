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

AVX2_32_SM_MT_datasets_gprof:
	gcc -O3 -pg $(src) SWAVX_TOP_datasets_MultiThread.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -lpthread -o SWAVX -D SUBMAT

AVX2_32_D_MT_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets_MultiThread.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -lpthread -o SWAVX -D DEBUG
AVX2_32_MT_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets_MultiThread.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -lpthread -o SWAVX


AVX2_WIN:
	x86_64-w64-mingw32-gcc -O3 $(src) SWAVX_TOP_datasets_MultiThread.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -lpthread -static -o SWAVX.exe
AVX2_WIN_BT:
	x86_64-w64-mingw32-gcc -O3 $(src) SWAVX_TOP_datasets_MultiThread.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -lpthread -static -o SWAVX.exe -D BT

AVX2_WIN_SM:
	x86_64-w64-mingw32-gcc -O3 $(src) SWAVX_TOP_datasets_MultiThread.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -lpthread -static -o SWAVX.exe -D SUBMAT
AVX2_WIN_SM_BT:
	x86_64-w64-mingw32-gcc -O3 $(src) SWAVX_TOP_datasets_MultiThread.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -lpthread -static -o SWAVX.exe -D SUBMAT -D BT

clean: 
	rm SWAVX




AVX2_16_D:
	gcc -O3 $(src) SWAVX_TOP.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -o SWAVX -D DEBUG -D L16 -D BT
AVX2_16_D_SM:
	gcc -O3 $(src) SWAVX_TOP.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -o SWAVX -D DEBUG -D L16 -D BT -D SUBMAT

AVX2_8_D:
	gcc -O3 $(src) SWAVX_TOP.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -o SWAVX -D DEBUG -D L8 -D BT
AVX2_8_D_SM:
	gcc -O3 $(src) SWAVX_TOP.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -o SWAVX -D DEBUG -D L8 -D BT -D SUBMAT


AVX2_16_SM_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -o SWAVX -D L16 -D SUBMAT
AVX2_16_SM_MT_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets_MultiThread.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -lpthread -o SWAVX -D L16 -D SUBMAT

AVX2_8_SM_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -o SWAVX -D L8 -D SUBMAT
AVX2_8_SM_MT_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets_MultiThread.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -lpthread -o SWAVX -D L8 -D SUBMAT

AVX2_16_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -o SWAVX -D L16
AVX2_16_MT_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets_MultiThread.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -lpthread -o SWAVX -D L16

AVX2_8_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -o SWAVX -D L8
AVX2_8_MT_datasets:
	gcc -O3 $(src) SWAVX_TOP_datasets_MultiThread.c SWAVX_256_Func_SeqToSeq_SubMat.c -mavx2 -lpthread -o SWAVX -D L8