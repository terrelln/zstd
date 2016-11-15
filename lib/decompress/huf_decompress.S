; ModuleID = '../lib/decompress/huf_decompress.c'
target datalayout = "e-m:o-i64:64-f80:128-n8:16:32:64-S128"
target triple = "x86_64-apple-macosx10.11.0"

%struct.algo_time_t = type { i32, i32 }
%struct.DTableDesc = type { i8, i8, i8, i8 }
%struct.HUF_DEltX2 = type { i8, i8 }
%struct.BIT_DStream_t = type { i64, i32, i8*, i8* }
%struct.sortedSymbol_t = type { i8, i8 }
%struct.HUF_DEltX4 = type { i16, i8, i8 }

@algoTime = internal unnamed_addr constant [16 x [3 x %struct.algo_time_t]] [[3 x %struct.algo_time_t] [%struct.algo_time_t zeroinitializer, %struct.algo_time_t { i32 1, i32 1 }, %struct.algo_time_t { i32 2, i32 2 }], [3 x %struct.algo_time_t] [%struct.algo_time_t zeroinitializer, %struct.algo_time_t { i32 1, i32 1 }, %struct.algo_time_t { i32 2, i32 2 }], [3 x %struct.algo_time_t] [%struct.algo_time_t { i32 38, i32 130 }, %struct.algo_time_t { i32 1313, i32 74 }, %struct.algo_time_t { i32 2151, i32 38 }], [3 x %struct.algo_time_t] [%struct.algo_time_t { i32 448, i32 128 }, %struct.algo_time_t { i32 1353, i32 74 }, %struct.algo_time_t { i32 2238, i32 41 }], [3 x %struct.algo_time_t] [%struct.algo_time_t { i32 556, i32 128 }, %struct.algo_time_t { i32 1353, i32 74 }, %struct.algo_time_t { i32 2238, i32 47 }], [3 x %struct.algo_time_t] [%struct.algo_time_t { i32 714, i32 128 }, %struct.algo_time_t { i32 1418, i32 74 }, %struct.algo_time_t { i32 2436, i32 53 }], [3 x %struct.algo_time_t] [%struct.algo_time_t { i32 883, i32 128 }, %struct.algo_time_t { i32 1437, i32 74 }, %struct.algo_time_t { i32 2464, i32 61 }], [3 x %struct.algo_time_t] [%struct.algo_time_t { i32 897, i32 128 }, %struct.algo_time_t { i32 1515, i32 75 }, %struct.algo_time_t { i32 2622, i32 68 }], [3 x %struct.algo_time_t] [%struct.algo_time_t { i32 926, i32 128 }, %struct.algo_time_t { i32 1613, i32 75 }, %struct.algo_time_t { i32 2730, i32 75 }], [3 x %struct.algo_time_t] [%struct.algo_time_t { i32 947, i32 128 }, %struct.algo_time_t { i32 1729, i32 77 }, %struct.algo_time_t { i32 3359, i32 77 }], [3 x %struct.algo_time_t] [%struct.algo_time_t { i32 1107, i32 128 }, %struct.algo_time_t { i32 2083, i32 81 }, %struct.algo_time_t { i32 4006, i32 84 }], [3 x %struct.algo_time_t] [%struct.algo_time_t { i32 1177, i32 128 }, %struct.algo_time_t { i32 2379, i32 87 }, %struct.algo_time_t { i32 4785, i32 88 }], [3 x %struct.algo_time_t] [%struct.algo_time_t { i32 1242, i32 128 }, %struct.algo_time_t { i32 2415, i32 93 }, %struct.algo_time_t { i32 5155, i32 84 }], [3 x %struct.algo_time_t] [%struct.algo_time_t { i32 1349, i32 128 }, %struct.algo_time_t { i32 2644, i32 106 }, %struct.algo_time_t { i32 5260, i32 106 }], [3 x %struct.algo_time_t] [%struct.algo_time_t { i32 1455, i32 128 }, %struct.algo_time_t { i32 2422, i32 124 }, %struct.algo_time_t { i32 4174, i32 124 }], [3 x %struct.algo_time_t] [%struct.algo_time_t { i32 722, i32 128 }, %struct.algo_time_t { i32 1891, i32 145 }, %struct.algo_time_t { i32 1936, i32 146 }]], align 16
@HUF_decompress.decompress = internal unnamed_addr constant [2 x i64 (i8*, i64, i8*, i64)*] [i64 (i8*, i64, i8*, i64)* @HUF_decompress4X2, i64 (i8*, i64, i8*, i64)* @HUF_decompress4X4], align 16

; Function Attrs: nounwind ssp uwtable
define i64 @HUF_readDTableX2(i32* %DTable, i8* %src, i64 %srcSize) #0 {
entry:
  %huffWeight = alloca [256 x i8], align 16
  %rankVal = alloca [17 x i32], align 16
  %tableLog = alloca i32, align 4
  %nbSymbols = alloca i32, align 4
  %dtd = alloca i32, align 4
  %tmpcast = bitcast i32* %dtd to %struct.DTableDesc*
  %0 = getelementptr inbounds [256 x i8], [256 x i8]* %huffWeight, i64 0, i64 0
  call void @llvm.lifetime.start(i64 256, i8* %0) #11
  %1 = bitcast [17 x i32]* %rankVal to i8*
  call void @llvm.lifetime.start(i64 68, i8* %1) #11
  %2 = bitcast i32* %tableLog to i8*
  call void @llvm.lifetime.start(i64 4, i8* %2) #11
  store i32 0, i32* %tableLog, align 4, !tbaa !2
  %3 = bitcast i32* %nbSymbols to i8*
  call void @llvm.lifetime.start(i64 4, i8* %3) #11
  store i32 0, i32* %nbSymbols, align 4, !tbaa !2
  %add.ptr = getelementptr inbounds i32, i32* %DTable, i64 1
  %4 = bitcast i32* %add.ptr to %struct.HUF_DEltX2*
  %arraydecay1 = getelementptr inbounds [17 x i32], [17 x i32]* %rankVal, i64 0, i64 0
  %call = call i64 @HUF_readStats(i8* %0, i64 256, i32* %arraydecay1, i32* nonnull %nbSymbols, i32* nonnull %tableLog, i8* %src, i64 %srcSize) #11
  %call2 = call i32 @HUF_isError(i64 %call) #11
  %tobool = icmp eq i32 %call2, 0
  br i1 %tobool, label %if.end, label %cleanup49

if.end:                                           ; preds = %entry
  %5 = bitcast i32* %dtd to i8*
  call void @llvm.lifetime.start(i64 4, i8* %5) #11
  %DTable.val = load i32, i32* %DTable, align 1
  store i32 %DTable.val, i32* %dtd, align 4
  %6 = load i32, i32* %tableLog, align 4, !tbaa !2
  %conv = and i32 %DTable.val, 255
  %add = add nuw nsw i32 %conv, 1
  %cmp = icmp ugt i32 %6, %add
  br i1 %cmp, label %cleanup, label %if.end6

if.end6:                                          ; preds = %if.end
  %tableType = getelementptr inbounds %struct.DTableDesc, %struct.DTableDesc* %tmpcast, i64 0, i32 1
  store i8 0, i8* %tableType, align 1, !tbaa !6
  %conv7 = trunc i32 %6 to i8
  %tableLog8 = getelementptr inbounds %struct.DTableDesc, %struct.DTableDesc* %tmpcast, i64 0, i32 2
  store i8 %conv7, i8* %tableLog8, align 2, !tbaa !8
  %7 = bitcast i32* %DTable to i8*
  %8 = call i64 @llvm.objectsize.i64.p0i8(i8* %7, i1 false)
  %call9 = call i8* @__memcpy_chk(i8* %7, i8* %5, i64 4, i64 %8) #11
  call void @llvm.lifetime.end(i64 4, i8* %5) #11
  %9 = load i32, i32* %tableLog, align 4, !tbaa !2
  %add1091 = add i32 %9, 1
  %cmp1192 = icmp ugt i32 %add1091, 1
  br i1 %cmp1192, label %for.body.lr.ph, label %for.cond17.preheader

for.body.lr.ph:                                   ; preds = %if.end6
  %10 = zext i32 %add1091 to i64
  %11 = add i32 %9, 1
  %12 = zext i32 %11 to i64
  %13 = add nuw nsw i64 %12, 3
  %14 = add nsw i64 %12, -2
  %xtraiter100 = and i64 %13, 3
  %lcmp.mod101 = icmp eq i64 %xtraiter100, 0
  br i1 %lcmp.mod101, label %for.body.lr.ph.split, label %for.body.prol.preheader

for.body.prol.preheader:                          ; preds = %for.body.lr.ph
  br label %for.body.prol

for.body.prol:                                    ; preds = %for.body.prol.preheader, %for.body.prol
  %indvars.iv97.prol = phi i64 [ %indvars.iv.next98.prol, %for.body.prol ], [ 1, %for.body.prol.preheader ]
  %nextRankStart.093.prol = phi i32 [ %add13.prol, %for.body.prol ], [ 0, %for.body.prol.preheader ]
  %prol.iter102 = phi i64 [ %prol.iter102.sub, %for.body.prol ], [ %xtraiter100, %for.body.prol.preheader ]
  %arrayidx.prol = getelementptr inbounds [17 x i32], [17 x i32]* %rankVal, i64 0, i64 %indvars.iv97.prol
  %15 = load i32, i32* %arrayidx.prol, align 4, !tbaa !2
  %16 = add i64 %indvars.iv97.prol, 4294967295
  %17 = trunc i64 %16 to i32
  %shl.prol = shl i32 %15, %17
  %add13.prol = add i32 %shl.prol, %nextRankStart.093.prol
  store i32 %nextRankStart.093.prol, i32* %arrayidx.prol, align 4, !tbaa !2
  %indvars.iv.next98.prol = add nuw nsw i64 %indvars.iv97.prol, 1
  %prol.iter102.sub = add i64 %prol.iter102, -1
  %prol.iter102.cmp = icmp eq i64 %prol.iter102.sub, 0
  br i1 %prol.iter102.cmp, label %for.body.lr.ph.split.loopexit, label %for.body.prol, !llvm.loop !9

for.body.lr.ph.split.loopexit:                    ; preds = %for.body.prol
  %indvars.iv.next98.prol.lcssa = phi i64 [ %indvars.iv.next98.prol, %for.body.prol ]
  %add13.prol.lcssa = phi i32 [ %add13.prol, %for.body.prol ]
  br label %for.body.lr.ph.split

for.body.lr.ph.split:                             ; preds = %for.body.lr.ph.split.loopexit, %for.body.lr.ph
  %indvars.iv97.unr = phi i64 [ 1, %for.body.lr.ph ], [ %indvars.iv.next98.prol.lcssa, %for.body.lr.ph.split.loopexit ]
  %nextRankStart.093.unr = phi i32 [ 0, %for.body.lr.ph ], [ %add13.prol.lcssa, %for.body.lr.ph.split.loopexit ]
  %18 = icmp ult i64 %14, 3
  br i1 %18, label %for.cond17.preheader.loopexit, label %for.body.lr.ph.split.split

for.body.lr.ph.split.split:                       ; preds = %for.body.lr.ph.split
  br label %for.body

cleanup:                                          ; preds = %if.end
  call void @llvm.lifetime.end(i64 4, i8* %5) #11
  br label %cleanup49

for.cond17.preheader.loopexit.unr-lcssa:          ; preds = %for.body
  br label %for.cond17.preheader.loopexit

for.cond17.preheader.loopexit:                    ; preds = %for.body.lr.ph.split, %for.cond17.preheader.loopexit.unr-lcssa
  br label %for.cond17.preheader

for.cond17.preheader:                             ; preds = %for.cond17.preheader.loopexit, %if.end6
  %19 = load i32, i32* %nbSymbols, align 4, !tbaa !2
  %cmp1889 = icmp eq i32 %19, 0
  br i1 %cmp1889, label %cleanup49, label %for.body20.lr.ph

for.body20.lr.ph:                                 ; preds = %for.cond17.preheader
  %20 = zext i32 %19 to i64
  br label %for.body20

for.body:                                         ; preds = %for.body, %for.body.lr.ph.split.split
  %indvars.iv97 = phi i64 [ %indvars.iv97.unr, %for.body.lr.ph.split.split ], [ %indvars.iv.next98.3, %for.body ]
  %nextRankStart.093 = phi i32 [ %nextRankStart.093.unr, %for.body.lr.ph.split.split ], [ %add13.3, %for.body ]
  %arrayidx = getelementptr inbounds [17 x i32], [17 x i32]* %rankVal, i64 0, i64 %indvars.iv97
  %21 = load i32, i32* %arrayidx, align 4, !tbaa !2
  %22 = add i64 %indvars.iv97, 4294967295
  %23 = trunc i64 %22 to i32
  %shl = shl i32 %21, %23
  %add13 = add i32 %shl, %nextRankStart.093
  store i32 %nextRankStart.093, i32* %arrayidx, align 4, !tbaa !2
  %indvars.iv.next98 = add nuw nsw i64 %indvars.iv97, 1
  %arrayidx.1 = getelementptr inbounds [17 x i32], [17 x i32]* %rankVal, i64 0, i64 %indvars.iv.next98
  %24 = load i32, i32* %arrayidx.1, align 4, !tbaa !2
  %25 = trunc i64 %indvars.iv97 to i32
  %shl.1 = shl i32 %24, %25
  %add13.1 = add i32 %shl.1, %add13
  store i32 %add13, i32* %arrayidx.1, align 4, !tbaa !2
  %indvars.iv.next98.1 = add nsw i64 %indvars.iv97, 2
  %arrayidx.2 = getelementptr inbounds [17 x i32], [17 x i32]* %rankVal, i64 0, i64 %indvars.iv.next98.1
  %26 = load i32, i32* %arrayidx.2, align 4, !tbaa !2
  %27 = add i64 %indvars.iv97, 1
  %28 = trunc i64 %27 to i32
  %shl.2 = shl i32 %26, %28
  %add13.2 = add i32 %shl.2, %add13.1
  store i32 %add13.1, i32* %arrayidx.2, align 4, !tbaa !2
  %indvars.iv.next98.2 = add nsw i64 %indvars.iv97, 3
  %arrayidx.3 = getelementptr inbounds [17 x i32], [17 x i32]* %rankVal, i64 0, i64 %indvars.iv.next98.2
  %29 = load i32, i32* %arrayidx.3, align 4, !tbaa !2
  %30 = add i64 %indvars.iv97, 2
  %31 = trunc i64 %30 to i32
  %shl.3 = shl i32 %29, %31
  %add13.3 = add i32 %shl.3, %add13.2
  store i32 %add13.2, i32* %arrayidx.3, align 4, !tbaa !2
  %indvars.iv.next98.3 = add nsw i64 %indvars.iv97, 4
  %cmp11.3 = icmp ult i64 %indvars.iv.next98.3, %10
  br i1 %cmp11.3, label %for.body, label %for.cond17.preheader.loopexit.unr-lcssa

for.body20:                                       ; preds = %for.body20.lr.ph, %for.end42
  %indvars.iv95 = phi i64 [ 0, %for.body20.lr.ph ], [ %indvars.iv.next96, %for.end42 ]
  %arrayidx22 = getelementptr inbounds [256 x i8], [256 x i8]* %huffWeight, i64 0, i64 %indvars.iv95
  %32 = load i8, i8* %arrayidx22, align 1, !tbaa !11
  %conv23 = zext i8 %32 to i32
  %shl24 = shl i32 1, %conv23
  %shr = ashr i32 %shl24, 1
  %conv25 = trunc i64 %indvars.iv95 to i8
  %sub27 = sub i32 %add1091, %conv23
  %conv28 = trunc i32 %sub27 to i8
  %idxprom29 = zext i8 %32 to i64
  %arrayidx30 = getelementptr inbounds [17 x i32], [17 x i32]* %rankVal, i64 0, i64 %idxprom29
  %33 = load i32, i32* %arrayidx30, align 4, !tbaa !2
  %add3486 = add i32 %33, %shr
  %cmp3587 = icmp ult i32 %33, %add3486
  br i1 %cmp3587, label %for.body37.lr.ph, label %for.end42

for.body37.lr.ph:                                 ; preds = %for.body20
  %34 = zext i32 %33 to i64
  %35 = zext i32 %add3486 to i64
  %36 = add i32 %33, %shr
  %37 = zext i32 %36 to i64
  %38 = zext i32 %33 to i64
  %39 = sub nsw i64 %37, %38
  %40 = add nsw i64 %37, -1
  %41 = sub nsw i64 %40, %38
  %xtraiter = and i64 %39, 7
  %lcmp.mod = icmp eq i64 %xtraiter, 0
  br i1 %lcmp.mod, label %for.body37.lr.ph.split, label %for.body37.prol.preheader

for.body37.prol.preheader:                        ; preds = %for.body37.lr.ph
  br label %for.body37.prol

for.body37.prol:                                  ; preds = %for.body37.prol, %for.body37.prol.preheader
  %indvars.iv.prol = phi i64 [ %indvars.iv.next.prol, %for.body37.prol ], [ %34, %for.body37.prol.preheader ]
  %prol.iter = phi i64 [ %prol.iter.sub, %for.body37.prol ], [ %xtraiter, %for.body37.prol.preheader ]
  %D.sroa.0.0..sroa_idx.prol = getelementptr inbounds %struct.HUF_DEltX2, %struct.HUF_DEltX2* %4, i64 %indvars.iv.prol, i32 0
  store i8 %conv25, i8* %D.sroa.0.0..sroa_idx.prol, align 1
  %D.sroa.4.0..sroa_idx.prol = getelementptr inbounds %struct.HUF_DEltX2, %struct.HUF_DEltX2* %4, i64 %indvars.iv.prol, i32 1
  store i8 %conv28, i8* %D.sroa.4.0..sroa_idx.prol, align 1
  %indvars.iv.next.prol = add nuw nsw i64 %indvars.iv.prol, 1
  %prol.iter.sub = add i64 %prol.iter, -1
  %prol.iter.cmp = icmp eq i64 %prol.iter.sub, 0
  br i1 %prol.iter.cmp, label %for.body37.lr.ph.split.loopexit, label %for.body37.prol, !llvm.loop !12

for.body37.lr.ph.split.loopexit:                  ; preds = %for.body37.prol
  %indvars.iv.next.prol.lcssa = phi i64 [ %indvars.iv.next.prol, %for.body37.prol ]
  br label %for.body37.lr.ph.split

for.body37.lr.ph.split:                           ; preds = %for.body37.lr.ph, %for.body37.lr.ph.split.loopexit
  %indvars.iv.unr = phi i64 [ %34, %for.body37.lr.ph ], [ %indvars.iv.next.prol.lcssa, %for.body37.lr.ph.split.loopexit ]
  %42 = icmp ult i64 %41, 7
  br i1 %42, label %for.end42.loopexit, label %for.body37.lr.ph.split.split

for.body37.lr.ph.split.split:                     ; preds = %for.body37.lr.ph.split
  br label %for.body37

for.body37:                                       ; preds = %for.body37, %for.body37.lr.ph.split.split
  %indvars.iv = phi i64 [ %indvars.iv.unr, %for.body37.lr.ph.split.split ], [ %indvars.iv.next.7, %for.body37 ]
  %D.sroa.0.0..sroa_idx = getelementptr inbounds %struct.HUF_DEltX2, %struct.HUF_DEltX2* %4, i64 %indvars.iv, i32 0
  store i8 %conv25, i8* %D.sroa.0.0..sroa_idx, align 1
  %D.sroa.4.0..sroa_idx = getelementptr inbounds %struct.HUF_DEltX2, %struct.HUF_DEltX2* %4, i64 %indvars.iv, i32 1
  store i8 %conv28, i8* %D.sroa.4.0..sroa_idx, align 1
  %indvars.iv.next = add nuw nsw i64 %indvars.iv, 1
  %D.sroa.0.0..sroa_idx.1 = getelementptr inbounds %struct.HUF_DEltX2, %struct.HUF_DEltX2* %4, i64 %indvars.iv.next, i32 0
  store i8 %conv25, i8* %D.sroa.0.0..sroa_idx.1, align 1
  %D.sroa.4.0..sroa_idx.1 = getelementptr inbounds %struct.HUF_DEltX2, %struct.HUF_DEltX2* %4, i64 %indvars.iv.next, i32 1
  store i8 %conv28, i8* %D.sroa.4.0..sroa_idx.1, align 1
  %indvars.iv.next.1 = add nsw i64 %indvars.iv, 2
  %D.sroa.0.0..sroa_idx.2 = getelementptr inbounds %struct.HUF_DEltX2, %struct.HUF_DEltX2* %4, i64 %indvars.iv.next.1, i32 0
  store i8 %conv25, i8* %D.sroa.0.0..sroa_idx.2, align 1
  %D.sroa.4.0..sroa_idx.2 = getelementptr inbounds %struct.HUF_DEltX2, %struct.HUF_DEltX2* %4, i64 %indvars.iv.next.1, i32 1
  store i8 %conv28, i8* %D.sroa.4.0..sroa_idx.2, align 1
  %indvars.iv.next.2 = add nsw i64 %indvars.iv, 3
  %D.sroa.0.0..sroa_idx.3 = getelementptr inbounds %struct.HUF_DEltX2, %struct.HUF_DEltX2* %4, i64 %indvars.iv.next.2, i32 0
  store i8 %conv25, i8* %D.sroa.0.0..sroa_idx.3, align 1
  %D.sroa.4.0..sroa_idx.3 = getelementptr inbounds %struct.HUF_DEltX2, %struct.HUF_DEltX2* %4, i64 %indvars.iv.next.2, i32 1
  store i8 %conv28, i8* %D.sroa.4.0..sroa_idx.3, align 1
  %indvars.iv.next.3 = add nsw i64 %indvars.iv, 4
  %D.sroa.0.0..sroa_idx.4 = getelementptr inbounds %struct.HUF_DEltX2, %struct.HUF_DEltX2* %4, i64 %indvars.iv.next.3, i32 0
  store i8 %conv25, i8* %D.sroa.0.0..sroa_idx.4, align 1
  %D.sroa.4.0..sroa_idx.4 = getelementptr inbounds %struct.HUF_DEltX2, %struct.HUF_DEltX2* %4, i64 %indvars.iv.next.3, i32 1
  store i8 %conv28, i8* %D.sroa.4.0..sroa_idx.4, align 1
  %indvars.iv.next.4 = add nsw i64 %indvars.iv, 5
  %D.sroa.0.0..sroa_idx.5 = getelementptr inbounds %struct.HUF_DEltX2, %struct.HUF_DEltX2* %4, i64 %indvars.iv.next.4, i32 0
  store i8 %conv25, i8* %D.sroa.0.0..sroa_idx.5, align 1
  %D.sroa.4.0..sroa_idx.5 = getelementptr inbounds %struct.HUF_DEltX2, %struct.HUF_DEltX2* %4, i64 %indvars.iv.next.4, i32 1
  store i8 %conv28, i8* %D.sroa.4.0..sroa_idx.5, align 1
  %indvars.iv.next.5 = add nsw i64 %indvars.iv, 6
  %D.sroa.0.0..sroa_idx.6 = getelementptr inbounds %struct.HUF_DEltX2, %struct.HUF_DEltX2* %4, i64 %indvars.iv.next.5, i32 0
  store i8 %conv25, i8* %D.sroa.0.0..sroa_idx.6, align 1
  %D.sroa.4.0..sroa_idx.6 = getelementptr inbounds %struct.HUF_DEltX2, %struct.HUF_DEltX2* %4, i64 %indvars.iv.next.5, i32 1
  store i8 %conv28, i8* %D.sroa.4.0..sroa_idx.6, align 1
  %indvars.iv.next.6 = add nsw i64 %indvars.iv, 7
  %D.sroa.0.0..sroa_idx.7 = getelementptr inbounds %struct.HUF_DEltX2, %struct.HUF_DEltX2* %4, i64 %indvars.iv.next.6, i32 0
  store i8 %conv25, i8* %D.sroa.0.0..sroa_idx.7, align 1
  %D.sroa.4.0..sroa_idx.7 = getelementptr inbounds %struct.HUF_DEltX2, %struct.HUF_DEltX2* %4, i64 %indvars.iv.next.6, i32 1
  store i8 %conv28, i8* %D.sroa.4.0..sroa_idx.7, align 1
  %indvars.iv.next.7 = add nsw i64 %indvars.iv, 8
  %cmp35.7 = icmp ult i64 %indvars.iv.next.7, %35
  br i1 %cmp35.7, label %for.body37, label %for.end42.loopexit.unr-lcssa

for.end42.loopexit.unr-lcssa:                     ; preds = %for.body37
  br label %for.end42.loopexit

for.end42.loopexit:                               ; preds = %for.body37.lr.ph.split, %for.end42.loopexit.unr-lcssa
  br label %for.end42

for.end42:                                        ; preds = %for.end42.loopexit, %for.body20
  store i32 %add3486, i32* %arrayidx30, align 4, !tbaa !2
  %indvars.iv.next96 = add nuw nsw i64 %indvars.iv95, 1
  %cmp18 = icmp ult i64 %indvars.iv.next96, %20
  br i1 %cmp18, label %for.body20, label %cleanup49.loopexit

cleanup49.loopexit:                               ; preds = %for.end42
  br label %cleanup49

cleanup49:                                        ; preds = %cleanup49.loopexit, %for.cond17.preheader, %cleanup, %entry
  %retval.1 = phi i64 [ -16, %cleanup ], [ %call, %entry ], [ %call, %for.cond17.preheader ], [ %call, %cleanup49.loopexit ]
  call void @llvm.lifetime.end(i64 4, i8* %3) #11
  call void @llvm.lifetime.end(i64 4, i8* %2) #11
  call void @llvm.lifetime.end(i64 68, i8* %1) #11
  call void @llvm.lifetime.end(i64 256, i8* %0) #11
  ret i64 %retval.1
}

; Function Attrs: argmemonly nounwind
declare void @llvm.lifetime.start(i64, i8* nocapture) #1

declare i64 @HUF_readStats(i8*, i64, i32*, i32*, i32*, i8*, i64) #2

declare i32 @HUF_isError(i64) #2

; Function Attrs: nounwind
declare i8* @__memcpy_chk(i8*, i8*, i64, i64) #3

; Function Attrs: nounwind readnone
declare i64 @llvm.objectsize.i64.p0i8(i8*, i1) #4

; Function Attrs: argmemonly nounwind
declare void @llvm.lifetime.end(i64, i8* nocapture) #1

; Function Attrs: argmemonly nounwind
declare void @llvm.memcpy.p0i8.p0i8.i64(i8* nocapture, i8* nocapture readonly, i64, i32, i1) #1

; Function Attrs: nounwind ssp uwtable
define i64 @HUF_decompress1X2_usingDTable(i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize, i32* nocapture readonly %DTable) #0 {
entry:
  %bitD.i = alloca %struct.BIT_DStream_t, align 8
  %DTable.val = load i32, i32* %DTable, align 1
  %conv = and i32 %DTable.val, 65280
  %cmp = icmp eq i32 %conv, 0
  br i1 %cmp, label %if.end, label %cleanup

if.end:                                           ; preds = %entry
  %0 = bitcast %struct.BIT_DStream_t* %bitD.i to i8*
  call void @llvm.lifetime.start(i64 32, i8* %0) #11
  %call2.i = call fastcc i64 @BIT_initDStream(%struct.BIT_DStream_t* nonnull %bitD.i, i8* %cSrc, i64 %cSrcSize) #11
  %call3.i = call i32 @HUF_isError(i64 %call2.i) #11
  %tobool.i = icmp eq i32 %call3.i, 0
  br i1 %tobool.i, label %cleanup.cont.i, label %HUF_decompress1X2_usingDTable_internal.exit

cleanup.cont.i:                                   ; preds = %if.end
  %dtd.sroa.3.0.extract.shift.i = lshr i32 %DTable.val, 16
  %conv.i = and i32 %dtd.sroa.3.0.extract.shift.i, 255
  %add.ptr1.i = getelementptr inbounds i32, i32* %DTable, i64 1
  %1 = bitcast i32* %add.ptr1.i to %struct.HUF_DEltX2*
  %add.ptr.i = getelementptr inbounds i8, i8* %dst, i64 %dstSize
  call fastcc void @HUF_decodeStreamX2(i8* %dst, %struct.BIT_DStream_t* nonnull %bitD.i, i8* %add.ptr.i, %struct.HUF_DEltX2* %1, i32 %conv.i) #11
  %bitD.idx.i = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD.i, i64 0, i32 1
  %bitD.idx.val.i = load i32, i32* %bitD.idx.i, align 8, !tbaa !13
  %bitD.idx28.i = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD.i, i64 0, i32 2
  %bitD.idx28.val.i = load i8*, i8** %bitD.idx28.i, align 8, !tbaa !17
  %bitD.idx29.i = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD.i, i64 0, i32 3
  %bitD.idx29.val.i = load i8*, i8** %bitD.idx29.i, align 8, !tbaa !18
  %cmp.i.i = icmp eq i8* %bitD.idx28.val.i, %bitD.idx29.val.i
  %cmp1.i.i = icmp eq i32 %bitD.idx.val.i, 64
  %tobool6.demorgan.i = and i1 %cmp1.i.i, %cmp.i.i
  %.dstSize.i = select i1 %tobool6.demorgan.i, i64 %dstSize, i64 -14
  br label %HUF_decompress1X2_usingDTable_internal.exit

HUF_decompress1X2_usingDTable_internal.exit:      ; preds = %if.end, %cleanup.cont.i
  %retval.1.i = phi i64 [ %call2.i, %if.end ], [ %.dstSize.i, %cleanup.cont.i ]
  call void @llvm.lifetime.end(i64 32, i8* %0) #11
  br label %cleanup

cleanup:                                          ; preds = %entry, %HUF_decompress1X2_usingDTable_internal.exit
  %retval.0 = phi i64 [ %retval.1.i, %HUF_decompress1X2_usingDTable_internal.exit ], [ -1, %entry ]
  ret i64 %retval.0
}

; Function Attrs: nounwind ssp uwtable
define i64 @HUF_decompress1X2_DCtx(i32* %DCtx, i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize) #0 {
entry:
  %bitD.i = alloca %struct.BIT_DStream_t, align 8
  %call = tail call i64 @HUF_readDTableX2(i32* %DCtx, i8* %cSrc, i64 %cSrcSize)
  %call1 = tail call i32 @HUF_isError(i64 %call) #11
  %tobool = icmp eq i32 %call1, 0
  br i1 %tobool, label %if.end, label %cleanup

if.end:                                           ; preds = %entry
  %cmp = icmp ult i64 %call, %cSrcSize
  br i1 %cmp, label %if.end3, label %cleanup

if.end3:                                          ; preds = %if.end
  %add.ptr = getelementptr inbounds i8, i8* %cSrc, i64 %call
  %sub = sub i64 %cSrcSize, %call
  %0 = bitcast %struct.BIT_DStream_t* %bitD.i to i8*
  call void @llvm.lifetime.start(i64 32, i8* %0) #11
  %DTable.val.i = load i32, i32* %DCtx, align 1
  %call2.i = call fastcc i64 @BIT_initDStream(%struct.BIT_DStream_t* nonnull %bitD.i, i8* %add.ptr, i64 %sub) #11
  %call3.i = call i32 @HUF_isError(i64 %call2.i) #11
  %tobool.i = icmp eq i32 %call3.i, 0
  br i1 %tobool.i, label %cleanup.cont.i, label %HUF_decompress1X2_usingDTable_internal.exit

cleanup.cont.i:                                   ; preds = %if.end3
  %dtd.sroa.3.0.extract.shift.i = lshr i32 %DTable.val.i, 16
  %conv.i = and i32 %dtd.sroa.3.0.extract.shift.i, 255
  %add.ptr1.i = getelementptr inbounds i32, i32* %DCtx, i64 1
  %1 = bitcast i32* %add.ptr1.i to %struct.HUF_DEltX2*
  %add.ptr.i = getelementptr inbounds i8, i8* %dst, i64 %dstSize
  call fastcc void @HUF_decodeStreamX2(i8* %dst, %struct.BIT_DStream_t* nonnull %bitD.i, i8* %add.ptr.i, %struct.HUF_DEltX2* %1, i32 %conv.i) #11
  %bitD.idx.i = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD.i, i64 0, i32 1
  %bitD.idx.val.i = load i32, i32* %bitD.idx.i, align 8, !tbaa !13
  %bitD.idx28.i = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD.i, i64 0, i32 2
  %bitD.idx28.val.i = load i8*, i8** %bitD.idx28.i, align 8, !tbaa !17
  %bitD.idx29.i = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD.i, i64 0, i32 3
  %bitD.idx29.val.i = load i8*, i8** %bitD.idx29.i, align 8, !tbaa !18
  %cmp.i.i = icmp eq i8* %bitD.idx28.val.i, %bitD.idx29.val.i
  %cmp1.i.i = icmp eq i32 %bitD.idx.val.i, 64
  %tobool6.demorgan.i = and i1 %cmp1.i.i, %cmp.i.i
  %.dstSize.i = select i1 %tobool6.demorgan.i, i64 %dstSize, i64 -14
  br label %HUF_decompress1X2_usingDTable_internal.exit

HUF_decompress1X2_usingDTable_internal.exit:      ; preds = %if.end3, %cleanup.cont.i
  %retval.1.i = phi i64 [ %call2.i, %if.end3 ], [ %.dstSize.i, %cleanup.cont.i ]
  call void @llvm.lifetime.end(i64 32, i8* %0) #11
  br label %cleanup

cleanup:                                          ; preds = %if.end, %entry, %HUF_decompress1X2_usingDTable_internal.exit
  %retval.0 = phi i64 [ %retval.1.i, %HUF_decompress1X2_usingDTable_internal.exit ], [ %call, %entry ], [ -13, %if.end ]
  ret i64 %retval.0
}

; Function Attrs: nounwind ssp uwtable
define i64 @HUF_decompress1X2(i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize) #0 {
entry:
  %DTable = alloca [2049 x i32], align 16
  %0 = bitcast [2049 x i32]* %DTable to i8*
  call void @llvm.lifetime.start(i64 8196, i8* %0) #11
  call void @llvm.memset.p0i8.i64(i8* %0, i8 0, i64 8196, i32 16, i1 false)
  %1 = getelementptr [2049 x i32], [2049 x i32]* %DTable, i64 0, i64 0
  store i32 184549387, i32* %1, align 16
  %call = call i64 @HUF_decompress1X2_DCtx(i32* %1, i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize)
  call void @llvm.lifetime.end(i64 8196, i8* %0) #11
  ret i64 %call
}

; Function Attrs: argmemonly nounwind
declare void @llvm.memset.p0i8.i64(i8* nocapture, i8, i64, i32, i1) #1

; Function Attrs: nounwind ssp uwtable
define i64 @HUF_decompress4X2_usingDTable(i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize, i32* nocapture readonly %DTable) #0 {
entry:
  %DTable.val = load i32, i32* %DTable, align 1
  %conv = and i32 %DTable.val, 65280
  %cmp = icmp eq i32 %conv, 0
  br i1 %cmp, label %if.end, label %cleanup

if.end:                                           ; preds = %entry
  %call2 = tail call fastcc i64 @HUF_decompress4X2_usingDTable_internal(i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize, i32* nonnull %DTable)
  br label %cleanup

cleanup:                                          ; preds = %entry, %if.end
  %retval.0 = phi i64 [ %call2, %if.end ], [ -1, %entry ]
  ret i64 %retval.0
}

; Function Attrs: nounwind ssp uwtable
define internal fastcc i64 @HUF_decompress4X2_usingDTable_internal(i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize, i32* nocapture readonly %DTable) unnamed_addr #0 {
entry:
  %bitD1 = alloca %struct.BIT_DStream_t, align 8
  %bitD2 = alloca %struct.BIT_DStream_t, align 8
  %bitD3 = alloca %struct.BIT_DStream_t, align 8
  %bitD4 = alloca %struct.BIT_DStream_t, align 8
  %cmp = icmp ult i64 %cSrcSize, 10
  br i1 %cmp, label %return, label %if.end

if.end:                                           ; preds = %entry
  %add.ptr = getelementptr inbounds i8, i8* %dst, i64 %dstSize
  %add.ptr1 = getelementptr inbounds i32, i32* %DTable, i64 1
  %0 = bitcast i32* %add.ptr1 to %struct.HUF_DEltX2*
  %1 = bitcast %struct.BIT_DStream_t* %bitD1 to i8*
  call void @llvm.lifetime.start(i64 32, i8* %1) #11
  %2 = bitcast %struct.BIT_DStream_t* %bitD2 to i8*
  call void @llvm.lifetime.start(i64 32, i8* %2) #11
  %3 = bitcast %struct.BIT_DStream_t* %bitD3 to i8*
  call void @llvm.lifetime.start(i64 32, i8* %3) #11
  %4 = bitcast %struct.BIT_DStream_t* %bitD4 to i8*
  call void @llvm.lifetime.start(i64 32, i8* %4) #11
  %val.0..sroa_cast1.i.i = bitcast i8* %cSrc to i16*
  %val.0.copyload.i.i = load i16, i16* %val.0..sroa_cast1.i.i, align 1
  %conv = zext i16 %val.0.copyload.i.i to i64
  %add.ptr2 = getelementptr inbounds i8, i8* %cSrc, i64 2
  %val.0..sroa_cast1.i.i362 = bitcast i8* %add.ptr2 to i16*
  %val.0.copyload.i.i363 = load i16, i16* %val.0..sroa_cast1.i.i362, align 1
  %conv4 = zext i16 %val.0.copyload.i.i363 to i64
  %add.ptr5 = getelementptr inbounds i8, i8* %cSrc, i64 4
  %val.0..sroa_cast1.i.i364 = bitcast i8* %add.ptr5 to i16*
  %val.0.copyload.i.i365 = load i16, i16* %val.0..sroa_cast1.i.i364, align 1
  %conv7 = zext i16 %val.0.copyload.i.i365 to i64
  %add.neg = add i64 %cSrcSize, -6
  %add8.neg = sub i64 %add.neg, %conv
  %add9.neg = sub i64 %add8.neg, %conv4
  %sub = sub i64 %add9.neg, %conv7
  %add.ptr10 = getelementptr inbounds i8, i8* %cSrc, i64 6
  %add.ptr11 = getelementptr inbounds i8, i8* %add.ptr10, i64 %conv
  %add.ptr12 = getelementptr inbounds i8, i8* %add.ptr11, i64 %conv4
  %add.ptr13 = getelementptr inbounds i8, i8* %add.ptr12, i64 %conv7
  %add14 = add i64 %dstSize, 3
  %div = lshr i64 %add14, 2
  %add.ptr15 = getelementptr inbounds i8, i8* %dst, i64 %div
  %add.ptr16 = getelementptr inbounds i8, i8* %add.ptr15, i64 %div
  %add.ptr17 = getelementptr inbounds i8, i8* %add.ptr16, i64 %div
  %DTable.val = load i32, i32* %DTable, align 1
  %dtd.sroa.3.0.extract.shift = lshr i32 %DTable.val, 16
  %conv19 = and i32 %dtd.sroa.3.0.extract.shift, 255
  %cmp20 = icmp ugt i64 %sub, %cSrcSize
  br i1 %cmp20, label %cleanup180, label %if.end23

if.end23:                                         ; preds = %if.end
  %call24 = call fastcc i64 @BIT_initDStream(%struct.BIT_DStream_t* nonnull %bitD1, i8* %add.ptr10, i64 %conv)
  %call25 = call i32 @HUF_isError(i64 %call24) #11
  %tobool = icmp eq i32 %call25, 0
  br i1 %tobool, label %cleanup.cont, label %cleanup180

cleanup.cont:                                     ; preds = %if.end23
  %call29 = call fastcc i64 @BIT_initDStream(%struct.BIT_DStream_t* nonnull %bitD2, i8* %add.ptr11, i64 %conv4)
  %call30 = call i32 @HUF_isError(i64 %call29) #11
  %tobool31 = icmp eq i32 %call30, 0
  br i1 %tobool31, label %cleanup.cont36, label %cleanup180

cleanup.cont36:                                   ; preds = %cleanup.cont
  %call38 = call fastcc i64 @BIT_initDStream(%struct.BIT_DStream_t* nonnull %bitD3, i8* %add.ptr12, i64 %conv7)
  %call39 = call i32 @HUF_isError(i64 %call38) #11
  %tobool40 = icmp eq i32 %call39, 0
  br i1 %tobool40, label %cleanup.cont45, label %cleanup180

cleanup.cont45:                                   ; preds = %cleanup.cont36
  %call47 = call fastcc i64 @BIT_initDStream(%struct.BIT_DStream_t* nonnull %bitD4, i8* %add.ptr13, i64 %sub)
  %call48 = call i32 @HUF_isError(i64 %call47) #11
  %tobool49 = icmp eq i32 %call48, 0
  br i1 %tobool49, label %cleanup.cont54, label %cleanup180

cleanup.cont54:                                   ; preds = %cleanup.cont45
  %call55 = call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* nonnull %bitD1)
  %call56 = call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* nonnull %bitD2)
  %or = or i32 %call56, %call55
  %call57 = call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* nonnull %bitD3)
  %or58 = or i32 %or, %call57
  %call59 = call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* nonnull %bitD4)
  %or60 = or i32 %or58, %call59
  %cmp61369 = icmp eq i32 %or60, 0
  %add.ptr63 = getelementptr inbounds i8, i8* %add.ptr, i64 -7
  %cmp64370 = icmp ult i8* %add.ptr17, %add.ptr63
  %or.cond371 = and i1 %cmp61369, %cmp64370
  br i1 %or.cond371, label %for.body.preheader, label %for.end

for.body.preheader:                               ; preds = %cleanup.cont54
  br label %for.body

for.body:                                         ; preds = %for.body.preheader, %for.body
  %op1.0375 = phi i8* [ %incdec.ptr141, %for.body ], [ %dst, %for.body.preheader ]
  %op2.0374 = phi i8* [ %incdec.ptr143, %for.body ], [ %add.ptr15, %for.body.preheader ]
  %op3.0373 = phi i8* [ %incdec.ptr145, %for.body ], [ %add.ptr16, %for.body.preheader ]
  %op4.0372 = phi i8* [ %incdec.ptr147, %for.body ], [ %add.ptr17, %for.body.preheader ]
  %call69 = call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* nonnull %bitD1, %struct.HUF_DEltX2* %0, i32 %conv19)
  %incdec.ptr = getelementptr inbounds i8, i8* %op1.0375, i64 1
  store i8 %call69, i8* %op1.0375, align 1, !tbaa !11
  %call74 = call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* nonnull %bitD2, %struct.HUF_DEltX2* %0, i32 %conv19)
  %incdec.ptr75 = getelementptr inbounds i8, i8* %op2.0374, i64 1
  store i8 %call74, i8* %op2.0374, align 1, !tbaa !11
  %call80 = call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* nonnull %bitD3, %struct.HUF_DEltX2* %0, i32 %conv19)
  %incdec.ptr81 = getelementptr inbounds i8, i8* %op3.0373, i64 1
  store i8 %call80, i8* %op3.0373, align 1, !tbaa !11
  %call86 = call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* nonnull %bitD4, %struct.HUF_DEltX2* %0, i32 %conv19)
  %incdec.ptr87 = getelementptr inbounds i8, i8* %op4.0372, i64 1
  store i8 %call86, i8* %op4.0372, align 1, !tbaa !11
  %call92 = call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* nonnull %bitD1, %struct.HUF_DEltX2* %0, i32 %conv19)
  %incdec.ptr93 = getelementptr inbounds i8, i8* %op1.0375, i64 2
  store i8 %call92, i8* %incdec.ptr, align 1, !tbaa !11
  %call99 = call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* nonnull %bitD2, %struct.HUF_DEltX2* %0, i32 %conv19)
  %incdec.ptr100 = getelementptr inbounds i8, i8* %op2.0374, i64 2
  store i8 %call99, i8* %incdec.ptr75, align 1, !tbaa !11
  %call106 = call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* nonnull %bitD3, %struct.HUF_DEltX2* %0, i32 %conv19)
  %incdec.ptr107 = getelementptr inbounds i8, i8* %op3.0373, i64 2
  store i8 %call106, i8* %incdec.ptr81, align 1, !tbaa !11
  %call113 = call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* nonnull %bitD4, %struct.HUF_DEltX2* %0, i32 %conv19)
  %incdec.ptr114 = getelementptr inbounds i8, i8* %op4.0372, i64 2
  store i8 %call113, i8* %incdec.ptr87, align 1, !tbaa !11
  %call119 = call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* nonnull %bitD1, %struct.HUF_DEltX2* %0, i32 %conv19)
  %incdec.ptr120 = getelementptr inbounds i8, i8* %op1.0375, i64 3
  store i8 %call119, i8* %incdec.ptr93, align 1, !tbaa !11
  %call125 = call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* nonnull %bitD2, %struct.HUF_DEltX2* %0, i32 %conv19)
  %incdec.ptr126 = getelementptr inbounds i8, i8* %op2.0374, i64 3
  store i8 %call125, i8* %incdec.ptr100, align 1, !tbaa !11
  %call131 = call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* nonnull %bitD3, %struct.HUF_DEltX2* %0, i32 %conv19)
  %incdec.ptr132 = getelementptr inbounds i8, i8* %op3.0373, i64 3
  store i8 %call131, i8* %incdec.ptr107, align 1, !tbaa !11
  %call137 = call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* nonnull %bitD4, %struct.HUF_DEltX2* %0, i32 %conv19)
  %incdec.ptr138 = getelementptr inbounds i8, i8* %op4.0372, i64 3
  store i8 %call137, i8* %incdec.ptr114, align 1, !tbaa !11
  %call140 = call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* nonnull %bitD1, %struct.HUF_DEltX2* %0, i32 %conv19)
  %incdec.ptr141 = getelementptr inbounds i8, i8* %op1.0375, i64 4
  store i8 %call140, i8* %incdec.ptr120, align 1, !tbaa !11
  %call142 = call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* nonnull %bitD2, %struct.HUF_DEltX2* %0, i32 %conv19)
  %incdec.ptr143 = getelementptr inbounds i8, i8* %op2.0374, i64 4
  store i8 %call142, i8* %incdec.ptr126, align 1, !tbaa !11
  %call144 = call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* nonnull %bitD3, %struct.HUF_DEltX2* %0, i32 %conv19)
  %incdec.ptr145 = getelementptr inbounds i8, i8* %op3.0373, i64 4
  store i8 %call144, i8* %incdec.ptr132, align 1, !tbaa !11
  %call146 = call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* nonnull %bitD4, %struct.HUF_DEltX2* %0, i32 %conv19)
  %incdec.ptr147 = getelementptr inbounds i8, i8* %op4.0372, i64 4
  store i8 %call146, i8* %incdec.ptr138, align 1, !tbaa !11
  %call148 = call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* nonnull %bitD1)
  %call149 = call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* nonnull %bitD2)
  %or150 = or i32 %call149, %call148
  %call151 = call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* nonnull %bitD3)
  %or152 = or i32 %or150, %call151
  %call153 = call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* nonnull %bitD4)
  %or154 = or i32 %or152, %call153
  %cmp61 = icmp eq i32 %or154, 0
  %cmp64 = icmp ult i8* %incdec.ptr147, %add.ptr63
  %or.cond = and i1 %cmp61, %cmp64
  br i1 %or.cond, label %for.body, label %for.end.loopexit

for.end.loopexit:                                 ; preds = %for.body
  %incdec.ptr147.lcssa = phi i8* [ %incdec.ptr147, %for.body ]
  %incdec.ptr145.lcssa = phi i8* [ %incdec.ptr145, %for.body ]
  %incdec.ptr143.lcssa = phi i8* [ %incdec.ptr143, %for.body ]
  %incdec.ptr141.lcssa = phi i8* [ %incdec.ptr141, %for.body ]
  br label %for.end

for.end:                                          ; preds = %for.end.loopexit, %cleanup.cont54
  %op1.0.lcssa = phi i8* [ %dst, %cleanup.cont54 ], [ %incdec.ptr141.lcssa, %for.end.loopexit ]
  %op2.0.lcssa = phi i8* [ %add.ptr15, %cleanup.cont54 ], [ %incdec.ptr143.lcssa, %for.end.loopexit ]
  %op3.0.lcssa = phi i8* [ %add.ptr16, %cleanup.cont54 ], [ %incdec.ptr145.lcssa, %for.end.loopexit ]
  %op4.0.lcssa = phi i8* [ %add.ptr17, %cleanup.cont54 ], [ %incdec.ptr147.lcssa, %for.end.loopexit ]
  %cmp155 = icmp ugt i8* %op1.0.lcssa, %add.ptr15
  %cmp159 = icmp ugt i8* %op2.0.lcssa, %add.ptr16
  %or.cond340 = or i1 %cmp159, %cmp155
  %cmp163 = icmp ugt i8* %op3.0.lcssa, %add.ptr17
  %or.cond341 = or i1 %cmp163, %or.cond340
  br i1 %or.cond341, label %cleanup180, label %if.end166

if.end166:                                        ; preds = %for.end
  call fastcc void @HUF_decodeStreamX2(i8* %op1.0.lcssa, %struct.BIT_DStream_t* nonnull %bitD1, i8* %add.ptr15, %struct.HUF_DEltX2* %0, i32 %conv19)
  call fastcc void @HUF_decodeStreamX2(i8* %op2.0.lcssa, %struct.BIT_DStream_t* nonnull %bitD2, i8* %add.ptr16, %struct.HUF_DEltX2* %0, i32 %conv19)
  call fastcc void @HUF_decodeStreamX2(i8* %op3.0.lcssa, %struct.BIT_DStream_t* nonnull %bitD3, i8* %add.ptr17, %struct.HUF_DEltX2* %0, i32 %conv19)
  call fastcc void @HUF_decodeStreamX2(i8* %op4.0.lcssa, %struct.BIT_DStream_t* nonnull %bitD4, i8* %add.ptr, %struct.HUF_DEltX2* %0, i32 %conv19)
  %bitD1.idx = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD1, i64 0, i32 1
  %bitD1.idx.val = load i32, i32* %bitD1.idx, align 8, !tbaa !13
  %bitD1.idx348 = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD1, i64 0, i32 2
  %bitD1.idx348.val = load i8*, i8** %bitD1.idx348, align 8, !tbaa !17
  %bitD1.idx349 = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD1, i64 0, i32 3
  %bitD1.idx349.val = load i8*, i8** %bitD1.idx349, align 8, !tbaa !18
  %cmp.i358 = icmp eq i8* %bitD1.idx348.val, %bitD1.idx349.val
  %cmp1.i359 = icmp eq i32 %bitD1.idx.val, 64
  %cmp1..i360 = and i1 %cmp1.i359, %cmp.i358
  %bitD2.idx = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD2, i64 0, i32 1
  %bitD2.idx.val = load i32, i32* %bitD2.idx, align 8, !tbaa !13
  %bitD2.idx346 = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD2, i64 0, i32 2
  %bitD2.idx346.val = load i8*, i8** %bitD2.idx346, align 8, !tbaa !17
  %bitD2.idx347 = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD2, i64 0, i32 3
  %bitD2.idx347.val = load i8*, i8** %bitD2.idx347, align 8, !tbaa !18
  %cmp.i354 = icmp eq i8* %bitD2.idx346.val, %bitD2.idx347.val
  %cmp1.i355 = icmp eq i32 %bitD2.idx.val, 64
  %cmp1..i356 = and i1 %cmp1.i355, %cmp.i354
  %and366 = and i1 %cmp1..i360, %cmp1..i356
  %bitD3.idx = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD3, i64 0, i32 1
  %bitD3.idx.val = load i32, i32* %bitD3.idx, align 8, !tbaa !13
  %bitD3.idx344 = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD3, i64 0, i32 2
  %bitD3.idx344.val = load i8*, i8** %bitD3.idx344, align 8, !tbaa !17
  %bitD3.idx345 = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD3, i64 0, i32 3
  %bitD3.idx345.val = load i8*, i8** %bitD3.idx345, align 8, !tbaa !18
  %cmp.i350 = icmp eq i8* %bitD3.idx344.val, %bitD3.idx345.val
  %cmp1.i351 = icmp eq i32 %bitD3.idx.val, 64
  %cmp1..i352 = and i1 %cmp1.i351, %cmp.i350
  %and174367 = and i1 %and366, %cmp1..i352
  %bitD4.idx = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD4, i64 0, i32 1
  %bitD4.idx.val = load i32, i32* %bitD4.idx, align 8, !tbaa !13
  %bitD4.idx342 = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD4, i64 0, i32 2
  %bitD4.idx342.val = load i8*, i8** %bitD4.idx342, align 8, !tbaa !17
  %bitD4.idx343 = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD4, i64 0, i32 3
  %bitD4.idx343.val = load i8*, i8** %bitD4.idx343, align 8, !tbaa !18
  %cmp.i = icmp eq i8* %bitD4.idx342.val, %bitD4.idx343.val
  %cmp1.i = icmp eq i32 %bitD4.idx.val, 64
  %cmp1..i = and i1 %cmp1.i, %cmp.i
  %and176368 = and i1 %and174367, %cmp1..i
  %.dstSize = select i1 %and176368, i64 %dstSize, i64 -14
  br label %cleanup180

cleanup180:                                       ; preds = %if.end166, %for.end, %if.end, %cleanup.cont45, %cleanup.cont36, %cleanup.cont, %if.end23
  %retval.4 = phi i64 [ %call47, %cleanup.cont45 ], [ %call38, %cleanup.cont36 ], [ %call29, %cleanup.cont ], [ %call24, %if.end23 ], [ -14, %if.end ], [ -14, %for.end ], [ %.dstSize, %if.end166 ]
  call void @llvm.lifetime.end(i64 32, i8* %4) #11
  call void @llvm.lifetime.end(i64 32, i8* %3) #11
  call void @llvm.lifetime.end(i64 32, i8* %2) #11
  call void @llvm.lifetime.end(i64 32, i8* %1) #11
  br label %return

return:                                           ; preds = %entry, %cleanup180
  %retval.5 = phi i64 [ %retval.4, %cleanup180 ], [ -14, %entry ]
  ret i64 %retval.5
}

; Function Attrs: nounwind ssp uwtable
define i64 @HUF_decompress4X2_DCtx(i32* %dctx, i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize) #0 {
entry:
  %call = tail call i64 @HUF_readDTableX2(i32* %dctx, i8* %cSrc, i64 %cSrcSize)
  %call1 = tail call i32 @HUF_isError(i64 %call) #11
  %tobool = icmp eq i32 %call1, 0
  br i1 %tobool, label %if.end, label %cleanup

if.end:                                           ; preds = %entry
  %cmp = icmp ult i64 %call, %cSrcSize
  br i1 %cmp, label %if.end3, label %cleanup

if.end3:                                          ; preds = %if.end
  %add.ptr = getelementptr inbounds i8, i8* %cSrc, i64 %call
  %sub = sub i64 %cSrcSize, %call
  %call4 = tail call fastcc i64 @HUF_decompress4X2_usingDTable_internal(i8* %dst, i64 %dstSize, i8* %add.ptr, i64 %sub, i32* %dctx)
  br label %cleanup

cleanup:                                          ; preds = %if.end, %entry, %if.end3
  %retval.0 = phi i64 [ %call4, %if.end3 ], [ %call, %entry ], [ -13, %if.end ]
  ret i64 %retval.0
}

; Function Attrs: nounwind ssp uwtable
define i64 @HUF_decompress4X2(i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize) #0 {
entry:
  %DTable = alloca [2049 x i32], align 16
  %0 = bitcast [2049 x i32]* %DTable to i8*
  call void @llvm.lifetime.start(i64 8196, i8* %0) #11
  call void @llvm.memset.p0i8.i64(i8* %0, i8 0, i64 8196, i32 16, i1 false)
  %1 = getelementptr [2049 x i32], [2049 x i32]* %DTable, i64 0, i64 0
  store i32 184549387, i32* %1, align 16
  %call.i = call i64 @HUF_readDTableX2(i32* %1, i8* %cSrc, i64 %cSrcSize) #11
  %call1.i = call i32 @HUF_isError(i64 %call.i) #11
  %tobool.i = icmp eq i32 %call1.i, 0
  br i1 %tobool.i, label %if.end.i, label %HUF_decompress4X2_DCtx.exit

if.end.i:                                         ; preds = %entry
  %cmp.i = icmp ult i64 %call.i, %cSrcSize
  br i1 %cmp.i, label %if.end3.i, label %HUF_decompress4X2_DCtx.exit

if.end3.i:                                        ; preds = %if.end.i
  %add.ptr.i = getelementptr inbounds i8, i8* %cSrc, i64 %call.i
  %sub.i = sub i64 %cSrcSize, %call.i
  %call4.i = call fastcc i64 @HUF_decompress4X2_usingDTable_internal(i8* %dst, i64 %dstSize, i8* %add.ptr.i, i64 %sub.i, i32* %1) #11
  br label %HUF_decompress4X2_DCtx.exit

HUF_decompress4X2_DCtx.exit:                      ; preds = %entry, %if.end.i, %if.end3.i
  %retval.0.i = phi i64 [ %call4.i, %if.end3.i ], [ %call.i, %entry ], [ -13, %if.end.i ]
  call void @llvm.lifetime.end(i64 8196, i8* %0) #11
  ret i64 %retval.0.i
}

; Function Attrs: nounwind ssp uwtable
define i64 @HUF_readDTableX4(i32* %DTable, i8* %src, i64 %srcSize) #0 {
entry:
  %rankVal.i.i = alloca [17 x i32], align 16
  %rankVal.i = alloca [17 x i32], align 16
  %weightList = alloca [256 x i8], align 16
  %sortedSymbol = alloca [256 x %struct.sortedSymbol_t], align 16
  %rankStats = alloca [17 x i32], align 16
  %rankStart0 = alloca [18 x i32], align 16
  %rankVal = alloca [16 x [17 x i32]], align 16
  %tableLog = alloca i32, align 4
  %nbSymbols = alloca i32, align 4
  %dtd = alloca i32, align 4
  %tmpcast = bitcast i32* %dtd to %struct.DTableDesc*
  %0 = getelementptr inbounds [256 x i8], [256 x i8]* %weightList, i64 0, i64 0
  call void @llvm.lifetime.start(i64 256, i8* %0) #11
  %1 = getelementptr inbounds [256 x %struct.sortedSymbol_t], [256 x %struct.sortedSymbol_t]* %sortedSymbol, i64 0, i64 0, i32 0
  call void @llvm.lifetime.start(i64 512, i8* %1) #11
  %2 = bitcast [17 x i32]* %rankStats to i8*
  call void @llvm.lifetime.start(i64 68, i8* %2) #11
  call void @llvm.memset.p0i8.i64(i8* %2, i8 0, i64 68, i32 16, i1 false)
  %3 = bitcast [18 x i32]* %rankStart0 to i8*
  call void @llvm.lifetime.start(i64 72, i8* %3) #11
  call void @llvm.memset.p0i8.i64(i8* %3, i8 0, i64 72, i32 16, i1 false)
  %add.ptr = getelementptr inbounds [18 x i32], [18 x i32]* %rankStart0, i64 0, i64 1
  %4 = bitcast [16 x [17 x i32]]* %rankVal to i8*
  call void @llvm.lifetime.start(i64 1088, i8* %4) #11
  %5 = bitcast i32* %tableLog to i8*
  call void @llvm.lifetime.start(i64 4, i8* %5) #11
  %6 = bitcast i32* %nbSymbols to i8*
  call void @llvm.lifetime.start(i64 4, i8* %6) #11
  %7 = bitcast i32* %dtd to i8*
  call void @llvm.lifetime.start(i64 4, i8* %7) #11
  %DTable.val = load i32, i32* %DTable, align 1
  store i32 %DTable.val, i32* %dtd, align 4
  %trunc = trunc i32 %DTable.val to i8
  %conv = and i32 %DTable.val, 255
  %add.ptr2 = getelementptr inbounds i32, i32* %DTable, i64 1
  %cmp = icmp ugt i8 %trunc, 16
  br i1 %cmp, label %cleanup, label %if.end

if.end:                                           ; preds = %entry
  %arraydecay5 = getelementptr inbounds [17 x i32], [17 x i32]* %rankStats, i64 0, i64 0
  %call6 = call i64 @HUF_readStats(i8* %0, i64 256, i32* %arraydecay5, i32* nonnull %nbSymbols, i32* nonnull %tableLog, i8* %src, i64 %srcSize) #11
  %call7 = call i32 @HUF_isError(i64 %call6) #11
  %tobool = icmp eq i32 %call7, 0
  br i1 %tobool, label %if.end9, label %cleanup

if.end9:                                          ; preds = %if.end
  %8 = load i32, i32* %tableLog, align 4, !tbaa !2
  %cmp10 = icmp ugt i32 %8, %conv
  br i1 %cmp10, label %cleanup, label %for.cond.preheader

for.cond.preheader:                               ; preds = %if.end9
  br label %for.cond

for.cond:                                         ; preds = %for.cond.preheader, %for.cond
  %maxW.0 = phi i32 [ %dec, %for.cond ], [ %8, %for.cond.preheader ]
  %idxprom = zext i32 %maxW.0 to i64
  %arrayidx = getelementptr inbounds [17 x i32], [17 x i32]* %rankStats, i64 0, i64 %idxprom
  %9 = load i32, i32* %arrayidx, align 4, !tbaa !2
  %cmp14 = icmp eq i32 %9, 0
  %dec = add i32 %maxW.0, -1
  br i1 %cmp14, label %for.cond, label %for.cond16.preheader

for.cond16.preheader:                             ; preds = %for.cond
  %maxW.0.lcssa = phi i32 [ %maxW.0, %for.cond ]
  %add = add i32 %maxW.0.lcssa, 1
  %cmp17195 = icmp ugt i32 %add, 1
  br i1 %cmp17195, label %for.body19.preheader, label %for.end26

for.body19.preheader:                             ; preds = %for.cond16.preheader
  %10 = add i32 %maxW.0.lcssa, -1
  %xtraiter290 = and i32 %maxW.0.lcssa, 3
  %lcmp.mod291 = icmp eq i32 %xtraiter290, 0
  br i1 %lcmp.mod291, label %for.body19.preheader.split, label %for.body19.prol.preheader

for.body19.prol.preheader:                        ; preds = %for.body19.preheader
  br label %for.body19.prol

for.body19.prol:                                  ; preds = %for.body19.prol.preheader, %for.body19.prol
  %indvars.iv216.prol = phi i64 [ %indvars.iv.next217.prol, %for.body19.prol ], [ 1, %for.body19.prol.preheader ]
  %nextRankStart.0196.prol = phi i32 [ %add22.prol, %for.body19.prol ], [ 0, %for.body19.prol.preheader ]
  %prol.iter292 = phi i32 [ %prol.iter292.sub, %for.body19.prol ], [ %xtraiter290, %for.body19.prol.preheader ]
  %arrayidx21.prol = getelementptr inbounds [17 x i32], [17 x i32]* %rankStats, i64 0, i64 %indvars.iv216.prol
  %11 = load i32, i32* %arrayidx21.prol, align 4, !tbaa !2
  %add22.prol = add i32 %11, %nextRankStart.0196.prol
  %arrayidx24.prol = getelementptr inbounds i32, i32* %add.ptr, i64 %indvars.iv216.prol
  store i32 %nextRankStart.0196.prol, i32* %arrayidx24.prol, align 4, !tbaa !2
  %indvars.iv.next217.prol = add nuw nsw i64 %indvars.iv216.prol, 1
  %prol.iter292.sub = add i32 %prol.iter292, -1
  %prol.iter292.cmp = icmp eq i32 %prol.iter292.sub, 0
  br i1 %prol.iter292.cmp, label %for.body19.preheader.split.loopexit, label %for.body19.prol, !llvm.loop !19

for.body19.preheader.split.loopexit:              ; preds = %for.body19.prol
  %indvars.iv.next217.prol.lcssa = phi i64 [ %indvars.iv.next217.prol, %for.body19.prol ]
  %add22.prol.lcssa = phi i32 [ %add22.prol, %for.body19.prol ]
  br label %for.body19.preheader.split

for.body19.preheader.split:                       ; preds = %for.body19.preheader.split.loopexit, %for.body19.preheader
  %add22.lcssa.unr = phi i32 [ undef, %for.body19.preheader ], [ %add22.prol.lcssa, %for.body19.preheader.split.loopexit ]
  %indvars.iv216.unr = phi i64 [ 1, %for.body19.preheader ], [ %indvars.iv.next217.prol.lcssa, %for.body19.preheader.split.loopexit ]
  %nextRankStart.0196.unr = phi i32 [ 0, %for.body19.preheader ], [ %add22.prol.lcssa, %for.body19.preheader.split.loopexit ]
  %12 = icmp ult i32 %10, 3
  br i1 %12, label %for.end26.loopexit, label %for.body19.preheader.split.split

for.body19.preheader.split.split:                 ; preds = %for.body19.preheader.split
  br label %for.body19

for.body19:                                       ; preds = %for.body19, %for.body19.preheader.split.split
  %indvars.iv216 = phi i64 [ %indvars.iv216.unr, %for.body19.preheader.split.split ], [ %indvars.iv.next217.3, %for.body19 ]
  %nextRankStart.0196 = phi i32 [ %nextRankStart.0196.unr, %for.body19.preheader.split.split ], [ %add22.3, %for.body19 ]
  %arrayidx21 = getelementptr inbounds [17 x i32], [17 x i32]* %rankStats, i64 0, i64 %indvars.iv216
  %13 = load i32, i32* %arrayidx21, align 4, !tbaa !2
  %add22 = add i32 %13, %nextRankStart.0196
  %arrayidx24 = getelementptr inbounds i32, i32* %add.ptr, i64 %indvars.iv216
  store i32 %nextRankStart.0196, i32* %arrayidx24, align 4, !tbaa !2
  %indvars.iv.next217 = add nuw nsw i64 %indvars.iv216, 1
  %arrayidx21.1 = getelementptr inbounds [17 x i32], [17 x i32]* %rankStats, i64 0, i64 %indvars.iv.next217
  %14 = load i32, i32* %arrayidx21.1, align 4, !tbaa !2
  %add22.1 = add i32 %14, %add22
  %arrayidx24.1 = getelementptr inbounds i32, i32* %add.ptr, i64 %indvars.iv.next217
  store i32 %add22, i32* %arrayidx24.1, align 4, !tbaa !2
  %indvars.iv.next217.1 = add nsw i64 %indvars.iv216, 2
  %arrayidx21.2 = getelementptr inbounds [17 x i32], [17 x i32]* %rankStats, i64 0, i64 %indvars.iv.next217.1
  %15 = load i32, i32* %arrayidx21.2, align 4, !tbaa !2
  %add22.2 = add i32 %15, %add22.1
  %arrayidx24.2 = getelementptr inbounds i32, i32* %add.ptr, i64 %indvars.iv.next217.1
  store i32 %add22.1, i32* %arrayidx24.2, align 4, !tbaa !2
  %indvars.iv.next217.2 = add nsw i64 %indvars.iv216, 3
  %arrayidx21.3 = getelementptr inbounds [17 x i32], [17 x i32]* %rankStats, i64 0, i64 %indvars.iv.next217.2
  %16 = load i32, i32* %arrayidx21.3, align 4, !tbaa !2
  %add22.3 = add i32 %16, %add22.2
  %arrayidx24.3 = getelementptr inbounds i32, i32* %add.ptr, i64 %indvars.iv.next217.2
  store i32 %add22.2, i32* %arrayidx24.3, align 4, !tbaa !2
  %indvars.iv.next217.3 = add nsw i64 %indvars.iv216, 4
  %lftr.wideiv218.3 = trunc i64 %indvars.iv.next217.3 to i32
  %exitcond219.3 = icmp eq i32 %lftr.wideiv218.3, %add
  br i1 %exitcond219.3, label %for.end26.loopexit.unr-lcssa, label %for.body19

for.end26.loopexit.unr-lcssa:                     ; preds = %for.body19
  %add22.3.lcssa = phi i32 [ %add22.3, %for.body19 ]
  br label %for.end26.loopexit

for.end26.loopexit:                               ; preds = %for.body19.preheader.split, %for.end26.loopexit.unr-lcssa
  %add22.lcssa = phi i32 [ %add22.lcssa.unr, %for.body19.preheader.split ], [ %add22.3.lcssa, %for.end26.loopexit.unr-lcssa ]
  br label %for.end26

for.end26:                                        ; preds = %for.end26.loopexit, %for.cond16.preheader
  %nextRankStart.0.lcssa = phi i32 [ 0, %for.cond16.preheader ], [ %add22.lcssa, %for.end26.loopexit ]
  store i32 %nextRankStart.0.lcssa, i32* %add.ptr, align 4, !tbaa !2
  %17 = load i32, i32* %nbSymbols, align 4, !tbaa !2
  %cmp29193 = icmp eq i32 %17, 0
  br i1 %cmp29193, label %for.end47, label %for.body31.lr.ph

for.body31.lr.ph:                                 ; preds = %for.end26
  %18 = zext i32 %17 to i64
  %xtraiter288 = and i32 %17, 1
  %lcmp.mod289 = icmp eq i32 %xtraiter288, 0
  br i1 %lcmp.mod289, label %for.body31.lr.ph.split, label %for.body31.prol

for.body31.prol:                                  ; preds = %for.body31.lr.ph
  %arrayidx34.prol = getelementptr inbounds [256 x i8], [256 x i8]* %weightList, i64 0, i64 0
  %19 = load i8, i8* %arrayidx34.prol, align 16, !tbaa !11
  %idxprom36.prol = zext i8 %19 to i64
  %arrayidx37.prol = getelementptr inbounds i32, i32* %add.ptr, i64 %idxprom36.prol
  %20 = load i32, i32* %arrayidx37.prol, align 4, !tbaa !2
  %inc38.prol = add i32 %20, 1
  store i32 %inc38.prol, i32* %arrayidx37.prol, align 4, !tbaa !2
  %idxprom40.prol = zext i32 %20 to i64
  %symbol.prol = getelementptr inbounds [256 x %struct.sortedSymbol_t], [256 x %struct.sortedSymbol_t]* %sortedSymbol, i64 0, i64 %idxprom40.prol, i32 0
  store i8 0, i8* %symbol.prol, align 2, !tbaa !20
  %weight.prol = getelementptr inbounds [256 x %struct.sortedSymbol_t], [256 x %struct.sortedSymbol_t]* %sortedSymbol, i64 0, i64 %idxprom40.prol, i32 1
  store i8 %19, i8* %weight.prol, align 1, !tbaa !22
  br label %for.body31.lr.ph.split

for.body31.lr.ph.split:                           ; preds = %for.body31.lr.ph, %for.body31.prol
  %indvars.iv214.unr = phi i64 [ 0, %for.body31.lr.ph ], [ 1, %for.body31.prol ]
  %21 = icmp eq i32 %17, 1
  br i1 %21, label %for.end47.loopexit, label %for.body31.lr.ph.split.split

for.body31.lr.ph.split.split:                     ; preds = %for.body31.lr.ph.split
  br label %for.body31

for.body31:                                       ; preds = %for.body31, %for.body31.lr.ph.split.split
  %indvars.iv214 = phi i64 [ %indvars.iv214.unr, %for.body31.lr.ph.split.split ], [ %indvars.iv.next215.1, %for.body31 ]
  %arrayidx34 = getelementptr inbounds [256 x i8], [256 x i8]* %weightList, i64 0, i64 %indvars.iv214
  %22 = load i8, i8* %arrayidx34, align 1, !tbaa !11
  %idxprom36 = zext i8 %22 to i64
  %arrayidx37 = getelementptr inbounds i32, i32* %add.ptr, i64 %idxprom36
  %23 = load i32, i32* %arrayidx37, align 4, !tbaa !2
  %inc38 = add i32 %23, 1
  store i32 %inc38, i32* %arrayidx37, align 4, !tbaa !2
  %conv39 = trunc i64 %indvars.iv214 to i8
  %idxprom40 = zext i32 %23 to i64
  %symbol = getelementptr inbounds [256 x %struct.sortedSymbol_t], [256 x %struct.sortedSymbol_t]* %sortedSymbol, i64 0, i64 %idxprom40, i32 0
  store i8 %conv39, i8* %symbol, align 2, !tbaa !20
  %weight = getelementptr inbounds [256 x %struct.sortedSymbol_t], [256 x %struct.sortedSymbol_t]* %sortedSymbol, i64 0, i64 %idxprom40, i32 1
  store i8 %22, i8* %weight, align 1, !tbaa !22
  %indvars.iv.next215 = add nuw nsw i64 %indvars.iv214, 1
  %arrayidx34.1 = getelementptr inbounds [256 x i8], [256 x i8]* %weightList, i64 0, i64 %indvars.iv.next215
  %24 = load i8, i8* %arrayidx34.1, align 1, !tbaa !11
  %idxprom36.1 = zext i8 %24 to i64
  %arrayidx37.1 = getelementptr inbounds i32, i32* %add.ptr, i64 %idxprom36.1
  %25 = load i32, i32* %arrayidx37.1, align 4, !tbaa !2
  %inc38.1 = add i32 %25, 1
  store i32 %inc38.1, i32* %arrayidx37.1, align 4, !tbaa !2
  %conv39.1 = trunc i64 %indvars.iv.next215 to i8
  %idxprom40.1 = zext i32 %25 to i64
  %symbol.1 = getelementptr inbounds [256 x %struct.sortedSymbol_t], [256 x %struct.sortedSymbol_t]* %sortedSymbol, i64 0, i64 %idxprom40.1, i32 0
  store i8 %conv39.1, i8* %symbol.1, align 2, !tbaa !20
  %weight.1 = getelementptr inbounds [256 x %struct.sortedSymbol_t], [256 x %struct.sortedSymbol_t]* %sortedSymbol, i64 0, i64 %idxprom40.1, i32 1
  store i8 %24, i8* %weight.1, align 1, !tbaa !22
  %indvars.iv.next215.1 = add nsw i64 %indvars.iv214, 2
  %cmp29.1 = icmp ult i64 %indvars.iv.next215.1, %18
  br i1 %cmp29.1, label %for.body31, label %for.end47.loopexit.unr-lcssa

for.end47.loopexit.unr-lcssa:                     ; preds = %for.body31
  br label %for.end47.loopexit

for.end47.loopexit:                               ; preds = %for.body31.lr.ph.split, %for.end47.loopexit.unr-lcssa
  br label %for.end47

for.end47:                                        ; preds = %for.end47.loopexit, %for.end26
  store i32 0, i32* %add.ptr, align 4, !tbaa !2
  br i1 %cmp17195, label %for.body57.lr.ph, label %for.end67

for.body57.lr.ph:                                 ; preds = %for.end47
  %sub = add nsw i32 %conv, -1
  %sub51 = sub i32 %sub, %8
  %26 = add i32 %maxW.0.lcssa, -1
  %xtraiter285 = and i32 %maxW.0.lcssa, 3
  %lcmp.mod286 = icmp eq i32 %xtraiter285, 0
  br i1 %lcmp.mod286, label %for.body57.lr.ph.split, label %for.body57.prol.preheader

for.body57.prol.preheader:                        ; preds = %for.body57.lr.ph
  br label %for.body57.prol

for.body57.prol:                                  ; preds = %for.body57.prol.preheader, %for.body57.prol
  %indvars.iv210.prol = phi i64 [ %indvars.iv.next211.prol, %for.body57.prol ], [ 1, %for.body57.prol.preheader ]
  %nextRankVal.0191.prol = phi i32 [ %add62.prol, %for.body57.prol ], [ 0, %for.body57.prol.preheader ]
  %prol.iter287 = phi i32 [ %prol.iter287.sub, %for.body57.prol ], [ %xtraiter285, %for.body57.prol.preheader ]
  %arrayidx60.prol = getelementptr inbounds [17 x i32], [17 x i32]* %rankStats, i64 0, i64 %indvars.iv210.prol
  %27 = load i32, i32* %arrayidx60.prol, align 4, !tbaa !2
  %28 = trunc i64 %indvars.iv210.prol to i32
  %add61.prol = add i32 %sub51, %28
  %shl.prol = shl i32 %27, %add61.prol
  %add62.prol = add i32 %shl.prol, %nextRankVal.0191.prol
  %arrayidx64.prol = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 0, i64 %indvars.iv210.prol
  store i32 %nextRankVal.0191.prol, i32* %arrayidx64.prol, align 4, !tbaa !2
  %indvars.iv.next211.prol = add nuw nsw i64 %indvars.iv210.prol, 1
  %prol.iter287.sub = add i32 %prol.iter287, -1
  %prol.iter287.cmp = icmp eq i32 %prol.iter287.sub, 0
  br i1 %prol.iter287.cmp, label %for.body57.lr.ph.split.loopexit, label %for.body57.prol, !llvm.loop !23

for.body57.lr.ph.split.loopexit:                  ; preds = %for.body57.prol
  %indvars.iv.next211.prol.lcssa = phi i64 [ %indvars.iv.next211.prol, %for.body57.prol ]
  %add62.prol.lcssa = phi i32 [ %add62.prol, %for.body57.prol ]
  br label %for.body57.lr.ph.split

for.body57.lr.ph.split:                           ; preds = %for.body57.lr.ph.split.loopexit, %for.body57.lr.ph
  %indvars.iv210.unr = phi i64 [ 1, %for.body57.lr.ph ], [ %indvars.iv.next211.prol.lcssa, %for.body57.lr.ph.split.loopexit ]
  %nextRankVal.0191.unr = phi i32 [ 0, %for.body57.lr.ph ], [ %add62.prol.lcssa, %for.body57.lr.ph.split.loopexit ]
  %29 = icmp ult i32 %26, 3
  br i1 %29, label %for.end67.loopexit, label %for.body57.lr.ph.split.split

for.body57.lr.ph.split.split:                     ; preds = %for.body57.lr.ph.split
  br label %for.body57

for.body57:                                       ; preds = %for.body57, %for.body57.lr.ph.split.split
  %indvars.iv210 = phi i64 [ %indvars.iv210.unr, %for.body57.lr.ph.split.split ], [ %indvars.iv.next211.3, %for.body57 ]
  %nextRankVal.0191 = phi i32 [ %nextRankVal.0191.unr, %for.body57.lr.ph.split.split ], [ %add62.3, %for.body57 ]
  %arrayidx60 = getelementptr inbounds [17 x i32], [17 x i32]* %rankStats, i64 0, i64 %indvars.iv210
  %30 = load i32, i32* %arrayidx60, align 4, !tbaa !2
  %31 = trunc i64 %indvars.iv210 to i32
  %add61 = add i32 %sub51, %31
  %shl = shl i32 %30, %add61
  %add62 = add i32 %shl, %nextRankVal.0191
  %arrayidx64 = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 0, i64 %indvars.iv210
  store i32 %nextRankVal.0191, i32* %arrayidx64, align 4, !tbaa !2
  %indvars.iv.next211 = add nuw nsw i64 %indvars.iv210, 1
  %arrayidx60.1 = getelementptr inbounds [17 x i32], [17 x i32]* %rankStats, i64 0, i64 %indvars.iv.next211
  %32 = load i32, i32* %arrayidx60.1, align 4, !tbaa !2
  %33 = trunc i64 %indvars.iv.next211 to i32
  %add61.1 = add i32 %sub51, %33
  %shl.1 = shl i32 %32, %add61.1
  %add62.1 = add i32 %shl.1, %add62
  %arrayidx64.1 = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 0, i64 %indvars.iv.next211
  store i32 %add62, i32* %arrayidx64.1, align 4, !tbaa !2
  %indvars.iv.next211.1 = add nsw i64 %indvars.iv210, 2
  %arrayidx60.2 = getelementptr inbounds [17 x i32], [17 x i32]* %rankStats, i64 0, i64 %indvars.iv.next211.1
  %34 = load i32, i32* %arrayidx60.2, align 4, !tbaa !2
  %35 = trunc i64 %indvars.iv.next211.1 to i32
  %add61.2 = add i32 %sub51, %35
  %shl.2 = shl i32 %34, %add61.2
  %add62.2 = add i32 %shl.2, %add62.1
  %arrayidx64.2 = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 0, i64 %indvars.iv.next211.1
  store i32 %add62.1, i32* %arrayidx64.2, align 4, !tbaa !2
  %indvars.iv.next211.2 = add nsw i64 %indvars.iv210, 3
  %arrayidx60.3 = getelementptr inbounds [17 x i32], [17 x i32]* %rankStats, i64 0, i64 %indvars.iv.next211.2
  %36 = load i32, i32* %arrayidx60.3, align 4, !tbaa !2
  %37 = trunc i64 %indvars.iv.next211.2 to i32
  %add61.3 = add i32 %sub51, %37
  %shl.3 = shl i32 %36, %add61.3
  %add62.3 = add i32 %shl.3, %add62.2
  %arrayidx64.3 = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 0, i64 %indvars.iv.next211.2
  store i32 %add62.2, i32* %arrayidx64.3, align 4, !tbaa !2
  %indvars.iv.next211.3 = add nsw i64 %indvars.iv210, 4
  %lftr.wideiv212.3 = trunc i64 %indvars.iv.next211.3 to i32
  %exitcond213.3 = icmp eq i32 %lftr.wideiv212.3, %add
  br i1 %exitcond213.3, label %for.end67.loopexit.unr-lcssa, label %for.body57

for.end67.loopexit.unr-lcssa:                     ; preds = %for.body57
  br label %for.end67.loopexit

for.end67.loopexit:                               ; preds = %for.body57.lr.ph.split, %for.end67.loopexit.unr-lcssa
  br label %for.end67

for.end67:                                        ; preds = %for.end67.loopexit, %for.end47
  %add68 = sub i32 1, %maxW.0.lcssa
  %sub69 = add i32 %add68, %8
  %sub71 = add nuw nsw i32 %conv, 1
  %add72 = sub i32 %sub71, %sub69
  %cmp73188.not = icmp uge i32 %sub69, %add72
  %cmp17195.not = xor i1 %cmp17195, true
  %brmerge = or i1 %cmp73188.not, %cmp17195.not
  br i1 %brmerge, label %for.end94, label %for.body75.us.preheader

for.body75.us.preheader:                          ; preds = %for.end67
  %38 = add i32 %8, 1
  %39 = sub i32 %38, %maxW.0.lcssa
  %40 = zext i32 %39 to i64
  %41 = add i32 %maxW.0.lcssa, %conv
  %42 = sub i32 %41, %8
  %43 = add i32 %maxW.0.lcssa, -1
  %44 = zext i32 %43 to i64
  %45 = add nuw nsw i64 %44, 1
  %46 = add i32 %8, 1
  %47 = sub i32 %46, %maxW.0.lcssa
  %48 = zext i32 %47 to i64
  %49 = zext i32 %43 to i64
  %50 = add nuw nsw i64 %49, 1
  %scevgep223 = getelementptr [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 0, i64 1
  %scevgep225 = getelementptr [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 0, i64 %50
  %51 = add i32 %maxW.0.lcssa, -1
  %52 = zext i32 %51 to i64
  %53 = add nuw nsw i64 %52, 1
  %54 = and i64 %53, 8589934588
  %55 = add nsw i64 %54, -4
  %56 = lshr exact i64 %55, 2
  %57 = add i32 %maxW.0.lcssa, 2
  %58 = add i32 %maxW.0.lcssa, 1
  %min.iters.check = icmp ult i64 %45, 4
  %n.vec = and i64 %45, 8589934588
  %cmp.zero = icmp eq i64 %n.vec, 0
  %ind.end = or i64 %n.vec, 1
  %59 = and i64 %56, 1
  %lcmp.mod281 = icmp eq i64 %59, 0
  %60 = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 0, i64 1
  %61 = bitcast i32* %60 to <2 x i32>*
  %62 = getelementptr [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 0, i64 3
  %63 = bitcast i32* %62 to <2 x i32>*
  %64 = icmp eq i64 %56, 0
  %cmp.n = icmp eq i64 %45, %n.vec
  br label %for.body75.us

for.body75.us:                                    ; preds = %for.cond80.for.end91_crit_edge.us, %for.body75.us.preheader
  %indvar = phi i64 [ %indvar.next, %for.cond80.for.end91_crit_edge.us ], [ 0, %for.body75.us.preheader ]
  %indvars.iv206 = phi i64 [ %indvars.iv.next207, %for.cond80.for.end91_crit_edge.us ], [ %40, %for.body75.us.preheader ]
  %consumed.0189.us = phi i32 [ %inc93.us, %for.cond80.for.end91_crit_edge.us ], [ %sub69, %for.body75.us.preheader ]
  %65 = add i64 %48, %indvar
  %scevgep = getelementptr [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 %65, i64 1
  %scevgep221 = getelementptr [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 %65, i64 %50
  %66 = trunc i64 %indvars.iv206 to i32
  br i1 %min.iters.check, label %for.body84.us.preheader, label %min.iters.checked

min.iters.checked:                                ; preds = %for.body75.us
  br i1 %cmp.zero, label %for.body84.us.preheader, label %vector.memcheck

vector.memcheck:                                  ; preds = %min.iters.checked
  %bound0 = icmp ule i32* %scevgep, %scevgep225
  %bound1 = icmp ule i32* %scevgep223, %scevgep221
  %memcheck.conflict = and i1 %bound0, %bound1
  br i1 %memcheck.conflict, label %for.body84.us.preheader, label %vector.ph

vector.ph:                                        ; preds = %vector.memcheck
  %broadcast.splatinsert229 = insertelement <2 x i32> undef, i32 %66, i32 0
  %broadcast.splat230 = shufflevector <2 x i32> %broadcast.splatinsert229, <2 x i32> undef, <2 x i32> zeroinitializer
  br i1 %lcmp.mod281, label %vector.body.prol, label %vector.ph.split

vector.body.prol:                                 ; preds = %vector.ph
  %wide.load.prol = load <2 x i32>, <2 x i32>* %61, align 4, !tbaa !2
  %wide.load228.prol = load <2 x i32>, <2 x i32>* %63, align 4, !tbaa !2
  %67 = lshr <2 x i32> %wide.load.prol, %broadcast.splat230
  %68 = lshr <2 x i32> %wide.load228.prol, %broadcast.splat230
  %69 = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 %indvars.iv206, i64 1
  %70 = bitcast i32* %69 to <2 x i32>*
  store <2 x i32> %67, <2 x i32>* %70, align 4, !tbaa !2
  %71 = getelementptr [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 %indvars.iv206, i64 3
  %72 = bitcast i32* %71 to <2 x i32>*
  store <2 x i32> %68, <2 x i32>* %72, align 4, !tbaa !2
  br label %vector.ph.split

vector.ph.split:                                  ; preds = %vector.body.prol, %vector.ph
  %index.unr = phi i64 [ 0, %vector.ph ], [ 4, %vector.body.prol ]
  br i1 %64, label %middle.block, label %vector.ph.split.split

vector.ph.split.split:                            ; preds = %vector.ph.split
  br label %vector.body

vector.body:                                      ; preds = %vector.body, %vector.ph.split.split
  %index = phi i64 [ %index.unr, %vector.ph.split.split ], [ %index.next.1, %vector.body ]
  %offset.idx = or i64 %index, 1
  %73 = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 0, i64 %offset.idx
  %74 = bitcast i32* %73 to <2 x i32>*
  %wide.load = load <2 x i32>, <2 x i32>* %74, align 4, !tbaa !2
  %75 = getelementptr i32, i32* %73, i64 2
  %76 = bitcast i32* %75 to <2 x i32>*
  %wide.load228 = load <2 x i32>, <2 x i32>* %76, align 4, !tbaa !2
  %77 = lshr <2 x i32> %wide.load, %broadcast.splat230
  %78 = lshr <2 x i32> %wide.load228, %broadcast.splat230
  %79 = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 %indvars.iv206, i64 %offset.idx
  %80 = bitcast i32* %79 to <2 x i32>*
  store <2 x i32> %77, <2 x i32>* %80, align 4, !tbaa !2
  %81 = getelementptr i32, i32* %79, i64 2
  %82 = bitcast i32* %81 to <2 x i32>*
  store <2 x i32> %78, <2 x i32>* %82, align 4, !tbaa !2
  %index.next = add i64 %index, 4
  %offset.idx.1 = or i64 %index.next, 1
  %83 = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 0, i64 %offset.idx.1
  %84 = bitcast i32* %83 to <2 x i32>*
  %wide.load.1 = load <2 x i32>, <2 x i32>* %84, align 4, !tbaa !2
  %85 = getelementptr i32, i32* %83, i64 2
  %86 = bitcast i32* %85 to <2 x i32>*
  %wide.load228.1 = load <2 x i32>, <2 x i32>* %86, align 4, !tbaa !2
  %87 = lshr <2 x i32> %wide.load.1, %broadcast.splat230
  %88 = lshr <2 x i32> %wide.load228.1, %broadcast.splat230
  %89 = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 %indvars.iv206, i64 %offset.idx.1
  %90 = bitcast i32* %89 to <2 x i32>*
  store <2 x i32> %87, <2 x i32>* %90, align 4, !tbaa !2
  %91 = getelementptr i32, i32* %89, i64 2
  %92 = bitcast i32* %91 to <2 x i32>*
  store <2 x i32> %88, <2 x i32>* %92, align 4, !tbaa !2
  %index.next.1 = add i64 %index, 8
  %93 = icmp eq i64 %index.next.1, %n.vec
  br i1 %93, label %middle.block.unr-lcssa, label %vector.body, !llvm.loop !24

middle.block.unr-lcssa:                           ; preds = %vector.body
  br label %middle.block

middle.block:                                     ; preds = %vector.ph.split, %middle.block.unr-lcssa
  br i1 %cmp.n, label %for.cond80.for.end91_crit_edge.us, label %for.body84.us.preheader

for.body84.us.preheader:                          ; preds = %middle.block, %vector.memcheck, %min.iters.checked, %for.body75.us
  %indvars.iv.ph = phi i64 [ 1, %vector.memcheck ], [ 1, %min.iters.checked ], [ 1, %for.body75.us ], [ %ind.end, %middle.block ]
  %94 = add nsw i64 %indvars.iv.ph, 1
  %95 = trunc i64 %94 to i32
  %96 = sub i32 %57, %95
  %97 = sub i32 %58, %95
  %xtraiter282 = and i32 %96, 3
  %lcmp.mod283 = icmp eq i32 %xtraiter282, 0
  br i1 %lcmp.mod283, label %for.body84.us.preheader.split, label %for.body84.us.prol.preheader

for.body84.us.prol.preheader:                     ; preds = %for.body84.us.preheader
  br label %for.body84.us.prol

for.body84.us.prol:                               ; preds = %for.body84.us.prol, %for.body84.us.prol.preheader
  %indvars.iv.prol = phi i64 [ %indvars.iv.next.prol, %for.body84.us.prol ], [ %indvars.iv.ph, %for.body84.us.prol.preheader ]
  %prol.iter284 = phi i32 [ %prol.iter284.sub, %for.body84.us.prol ], [ %xtraiter282, %for.body84.us.prol.preheader ]
  %arrayidx86.us.prol = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 0, i64 %indvars.iv.prol
  %98 = load i32, i32* %arrayidx86.us.prol, align 4, !tbaa !2
  %shr.us.prol = lshr i32 %98, %66
  %arrayidx88.us.prol = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 %indvars.iv206, i64 %indvars.iv.prol
  store i32 %shr.us.prol, i32* %arrayidx88.us.prol, align 4, !tbaa !2
  %indvars.iv.next.prol = add nuw nsw i64 %indvars.iv.prol, 1
  %prol.iter284.sub = add i32 %prol.iter284, -1
  %prol.iter284.cmp = icmp eq i32 %prol.iter284.sub, 0
  br i1 %prol.iter284.cmp, label %for.body84.us.preheader.split.loopexit, label %for.body84.us.prol, !llvm.loop !27

for.body84.us.preheader.split.loopexit:           ; preds = %for.body84.us.prol
  %indvars.iv.next.prol.lcssa = phi i64 [ %indvars.iv.next.prol, %for.body84.us.prol ]
  br label %for.body84.us.preheader.split

for.body84.us.preheader.split:                    ; preds = %for.body84.us.preheader, %for.body84.us.preheader.split.loopexit
  %indvars.iv.unr = phi i64 [ %indvars.iv.ph, %for.body84.us.preheader ], [ %indvars.iv.next.prol.lcssa, %for.body84.us.preheader.split.loopexit ]
  %99 = icmp ult i32 %97, 3
  br i1 %99, label %for.cond80.for.end91_crit_edge.us.loopexit, label %for.body84.us.preheader.split.split

for.body84.us.preheader.split.split:              ; preds = %for.body84.us.preheader.split
  br label %for.body84.us

for.body84.us:                                    ; preds = %for.body84.us, %for.body84.us.preheader.split.split
  %indvars.iv = phi i64 [ %indvars.iv.unr, %for.body84.us.preheader.split.split ], [ %indvars.iv.next.3, %for.body84.us ]
  %arrayidx86.us = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 0, i64 %indvars.iv
  %100 = load i32, i32* %arrayidx86.us, align 4, !tbaa !2
  %shr.us = lshr i32 %100, %66
  %arrayidx88.us = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 %indvars.iv206, i64 %indvars.iv
  store i32 %shr.us, i32* %arrayidx88.us, align 4, !tbaa !2
  %indvars.iv.next = add nuw nsw i64 %indvars.iv, 1
  %arrayidx86.us.1 = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 0, i64 %indvars.iv.next
  %101 = load i32, i32* %arrayidx86.us.1, align 4, !tbaa !2
  %shr.us.1 = lshr i32 %101, %66
  %arrayidx88.us.1 = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 %indvars.iv206, i64 %indvars.iv.next
  store i32 %shr.us.1, i32* %arrayidx88.us.1, align 4, !tbaa !2
  %indvars.iv.next.1 = add nsw i64 %indvars.iv, 2
  %arrayidx86.us.2 = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 0, i64 %indvars.iv.next.1
  %102 = load i32, i32* %arrayidx86.us.2, align 4, !tbaa !2
  %shr.us.2 = lshr i32 %102, %66
  %arrayidx88.us.2 = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 %indvars.iv206, i64 %indvars.iv.next.1
  store i32 %shr.us.2, i32* %arrayidx88.us.2, align 4, !tbaa !2
  %indvars.iv.next.2 = add nsw i64 %indvars.iv, 3
  %arrayidx86.us.3 = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 0, i64 %indvars.iv.next.2
  %103 = load i32, i32* %arrayidx86.us.3, align 4, !tbaa !2
  %shr.us.3 = lshr i32 %103, %66
  %arrayidx88.us.3 = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 %indvars.iv206, i64 %indvars.iv.next.2
  store i32 %shr.us.3, i32* %arrayidx88.us.3, align 4, !tbaa !2
  %indvars.iv.next.3 = add nsw i64 %indvars.iv, 4
  %lftr.wideiv204.3 = trunc i64 %indvars.iv.next.3 to i32
  %exitcond205.3 = icmp eq i32 %lftr.wideiv204.3, %add
  br i1 %exitcond205.3, label %for.cond80.for.end91_crit_edge.us.loopexit.unr-lcssa, label %for.body84.us, !llvm.loop !28

for.cond80.for.end91_crit_edge.us.loopexit.unr-lcssa: ; preds = %for.body84.us
  br label %for.cond80.for.end91_crit_edge.us.loopexit

for.cond80.for.end91_crit_edge.us.loopexit:       ; preds = %for.body84.us.preheader.split, %for.cond80.for.end91_crit_edge.us.loopexit.unr-lcssa
  br label %for.cond80.for.end91_crit_edge.us

for.cond80.for.end91_crit_edge.us:                ; preds = %for.cond80.for.end91_crit_edge.us.loopexit, %middle.block
  %inc93.us = add i32 %consumed.0189.us, 1
  %indvars.iv.next207 = add nuw nsw i64 %indvars.iv206, 1
  %exitcond208 = icmp eq i32 %inc93.us, %42
  %indvar.next = add i64 %indvar, 1
  br i1 %exitcond208, label %for.end94.loopexit, label %for.body75.us

for.end94.loopexit:                               ; preds = %for.cond80.for.end91_crit_edge.us
  br label %for.end94

for.end94:                                        ; preds = %for.end94.loopexit, %for.end67
  %add98 = add i32 %8, 1
  %104 = bitcast [17 x i32]* %rankVal.i to i8*
  call void @llvm.lifetime.start(i64 68, i8* %104) #11
  %sub.i = sub i32 %add98, %conv
  %sub1.i = sub i32 %add98, %maxW.0.lcssa
  call void @llvm.memcpy.p0i8.p0i8.i64(i8* %104, i8* %4, i64 68, i32 16, i1 false) #11
  %cmp103.i = icmp eq i32 %nextRankStart.0.lcssa, 0
  br i1 %cmp103.i, label %HUF_fillDTableX4.exit, label %for.body.lr.ph.i

for.body.lr.ph.i:                                 ; preds = %for.end94
  %105 = bitcast [17 x i32]* %rankVal.i.i to i8*
  br label %for.body.i

for.body.i:                                       ; preds = %if.end36.i, %for.body.lr.ph.i
  %indvars.iv109.i = phi i64 [ 0, %for.body.lr.ph.i ], [ %indvars.iv.next110.i, %if.end36.i ]
  %symbol2.i = getelementptr inbounds [256 x %struct.sortedSymbol_t], [256 x %struct.sortedSymbol_t]* %sortedSymbol, i64 0, i64 %indvars.iv109.i, i32 0
  %106 = load i8, i8* %symbol2.i, align 2, !tbaa !20
  %conv.i = zext i8 %106 to i32
  %weight5.i = getelementptr inbounds [256 x %struct.sortedSymbol_t], [256 x %struct.sortedSymbol_t]* %sortedSymbol, i64 0, i64 %indvars.iv109.i, i32 1
  %107 = load i8, i8* %weight5.i, align 1, !tbaa !22
  %conv6.i = zext i8 %107 to i32
  %sub7.i = sub i32 %add98, %conv6.i
  %idxprom8.i = zext i8 %107 to i64
  %arrayidx9.i = getelementptr inbounds [17 x i32], [17 x i32]* %rankVal.i, i64 0, i64 %idxprom8.i
  %108 = load i32, i32* %arrayidx9.i, align 4, !tbaa !2
  %sub10.i = sub i32 %conv, %sub7.i
  %shl.i = shl i32 1, %sub10.i
  %cmp12.i = icmp ult i32 %sub10.i, %sub1.i
  br i1 %cmp12.i, label %if.else.i, label %if.then.i

if.then.i:                                        ; preds = %for.body.i
  %add.i = add i32 %sub7.i, %sub.i
  %cmp14.i = icmp slt i32 %add.i, 1
  %.add.i = select i1 %cmp14.i, i32 1, i32 %add.i
  %idxprom17.i = sext i32 %.add.i to i64
  %arrayidx18.i = getelementptr inbounds [18 x i32], [18 x i32]* %rankStart0, i64 0, i64 %idxprom17.i
  %109 = load i32, i32* %arrayidx18.i, align 4, !tbaa !2
  %idx.ext.i = zext i32 %108 to i64
  %add.ptr.i182 = getelementptr inbounds i32, i32* %add.ptr2, i64 %idx.ext.i
  %idxprom20.i = zext i32 %sub7.i to i64
  %arraydecay22.i = getelementptr inbounds [16 x [17 x i32]], [16 x [17 x i32]]* %rankVal, i64 0, i64 %idxprom20.i, i64 0
  %idx.ext23.i = zext i32 %109 to i64
  %add.ptr24.i = getelementptr inbounds [256 x %struct.sortedSymbol_t], [256 x %struct.sortedSymbol_t]* %sortedSymbol, i64 0, i64 %idx.ext23.i
  call void @llvm.lifetime.start(i64 68, i8* %105) #11
  %110 = bitcast i32* %arraydecay22.i to i8*
  call void @llvm.memcpy.p0i8.p0i8.i64(i8* %105, i8* %110, i64 68, i32 4, i1 false) #11
  %cmp.i.i = icmp sgt i32 %.add.i, 1
  br i1 %cmp.i.i, label %if.then.i.i, label %for.cond5.preheader.i.i

for.cond5.preheader.i.i.loopexit:                 ; preds = %for.body.i.i
  br label %for.cond5.preheader.i.i

for.cond5.preheader.i.i:                          ; preds = %for.cond5.preheader.i.i.loopexit, %middle.block257, %if.then.i.i, %if.then.i
  %cmp671.i.i = icmp eq i32 %nextRankStart.0.lcssa, %109
  br i1 %cmp671.i.i, label %HUF_fillDTableX4Level2.exit.i, label %for.body8.i.preheader.i

for.body8.i.preheader.i:                          ; preds = %for.cond5.preheader.i.i
  %111 = sub i32 %nextRankStart.0.lcssa, %109
  br label %for.body8.i.i

if.then.i.i:                                      ; preds = %if.then.i
  %arrayidx.i.i = getelementptr inbounds [17 x i32], [17 x i32]* %rankVal.i.i, i64 0, i64 %idxprom17.i
  %112 = load i32, i32* %arrayidx.i.i, align 4, !tbaa !2
  %cmp173.i.i = icmp eq i32 %112, 0
  br i1 %cmp173.i.i, label %for.cond5.preheader.i.i, label %for.body.lr.ph.i.i

for.body.lr.ph.i.i:                               ; preds = %if.then.i.i
  %DElt.i.sroa.6.0.insert.ext92.i = shl i32 %sub7.i, 16
  %DElt.i.sroa.6.0.insert.shift93.i = and i32 %DElt.i.sroa.6.0.insert.ext92.i, 16711680
  %DElt.i.sroa.6.0.insert.insert95.i = or i32 %conv.i, %DElt.i.sroa.6.0.insert.shift93.i
  %DElt.i.sroa.0.0.insert.insert90.i = or i32 %DElt.i.sroa.6.0.insert.insert95.i, 16777216
  %113 = add i32 %112, -1
  %114 = zext i32 %113 to i64
  %115 = add nuw nsw i64 %114, 1
  %min.iters.check259 = icmp ult i64 %115, 8
  br i1 %min.iters.check259, label %for.body.i.i.preheader, label %min.iters.checked260

min.iters.checked260:                             ; preds = %for.body.lr.ph.i.i
  %n.vec262 = and i64 %115, 8589934584
  %cmp.zero263 = icmp eq i64 %n.vec262, 0
  br i1 %cmp.zero263, label %for.body.i.i.preheader, label %vector.ph264

vector.ph264:                                     ; preds = %min.iters.checked260
  %broadcast.splatinsert273 = insertelement <4 x i32> undef, i32 %DElt.i.sroa.0.0.insert.insert90.i, i32 0
  %broadcast.splat274 = shufflevector <4 x i32> %broadcast.splatinsert273, <4 x i32> undef, <4 x i32> zeroinitializer
  %116 = add i32 %112, -1
  %117 = zext i32 %116 to i64
  %118 = add nuw nsw i64 %117, 1
  %119 = and i64 %118, 8589934584
  %120 = add nsw i64 %119, -8
  %121 = lshr exact i64 %120, 3
  %122 = add nuw nsw i64 %121, 1
  %xtraiter = and i64 %122, 7
  %lcmp.mod = icmp eq i64 %xtraiter, 0
  br i1 %lcmp.mod, label %vector.ph264.split, label %vector.body256.prol.preheader

vector.body256.prol.preheader:                    ; preds = %vector.ph264
  br label %vector.body256.prol

vector.body256.prol:                              ; preds = %vector.body256.prol, %vector.body256.prol.preheader
  %index265.prol = phi i64 [ %index.next266.prol, %vector.body256.prol ], [ 0, %vector.body256.prol.preheader ]
  %prol.iter = phi i64 [ %prol.iter.sub, %vector.body256.prol ], [ %xtraiter, %vector.body256.prol.preheader ]
  %123 = getelementptr inbounds i32, i32* %add.ptr.i182, i64 %index265.prol
  %124 = bitcast i32* %123 to <4 x i32>*
  store <4 x i32> %broadcast.splat274, <4 x i32>* %124, align 2
  %125 = getelementptr i32, i32* %123, i64 4
  %126 = bitcast i32* %125 to <4 x i32>*
  store <4 x i32> %broadcast.splat274, <4 x i32>* %126, align 2
  %index.next266.prol = add i64 %index265.prol, 8
  %prol.iter.sub = add i64 %prol.iter, -1
  %prol.iter.cmp = icmp eq i64 %prol.iter.sub, 0
  br i1 %prol.iter.cmp, label %vector.ph264.split.loopexit, label %vector.body256.prol, !llvm.loop !29

vector.ph264.split.loopexit:                      ; preds = %vector.body256.prol
  %index.next266.prol.lcssa = phi i64 [ %index.next266.prol, %vector.body256.prol ]
  br label %vector.ph264.split

vector.ph264.split:                               ; preds = %vector.ph264, %vector.ph264.split.loopexit
  %index265.unr = phi i64 [ 0, %vector.ph264 ], [ %index.next266.prol.lcssa, %vector.ph264.split.loopexit ]
  %127 = icmp ult i64 %120, 56
  br i1 %127, label %middle.block257, label %vector.ph264.split.split

vector.ph264.split.split:                         ; preds = %vector.ph264.split
  br label %vector.body256

vector.body256:                                   ; preds = %vector.body256, %vector.ph264.split.split
  %index265 = phi i64 [ %index265.unr, %vector.ph264.split.split ], [ %index.next266.7, %vector.body256 ]
  %128 = getelementptr inbounds i32, i32* %add.ptr.i182, i64 %index265
  %129 = bitcast i32* %128 to <4 x i32>*
  store <4 x i32> %broadcast.splat274, <4 x i32>* %129, align 2
  %130 = getelementptr i32, i32* %128, i64 4
  %131 = bitcast i32* %130 to <4 x i32>*
  store <4 x i32> %broadcast.splat274, <4 x i32>* %131, align 2
  %index.next266 = add i64 %index265, 8
  %132 = getelementptr inbounds i32, i32* %add.ptr.i182, i64 %index.next266
  %133 = bitcast i32* %132 to <4 x i32>*
  store <4 x i32> %broadcast.splat274, <4 x i32>* %133, align 2
  %134 = getelementptr i32, i32* %132, i64 4
  %135 = bitcast i32* %134 to <4 x i32>*
  store <4 x i32> %broadcast.splat274, <4 x i32>* %135, align 2
  %index.next266.1 = add i64 %index265, 16
  %136 = getelementptr inbounds i32, i32* %add.ptr.i182, i64 %index.next266.1
  %137 = bitcast i32* %136 to <4 x i32>*
  store <4 x i32> %broadcast.splat274, <4 x i32>* %137, align 2
  %138 = getelementptr i32, i32* %136, i64 4
  %139 = bitcast i32* %138 to <4 x i32>*
  store <4 x i32> %broadcast.splat274, <4 x i32>* %139, align 2
  %index.next266.2 = add i64 %index265, 24
  %140 = getelementptr inbounds i32, i32* %add.ptr.i182, i64 %index.next266.2
  %141 = bitcast i32* %140 to <4 x i32>*
  store <4 x i32> %broadcast.splat274, <4 x i32>* %141, align 2
  %142 = getelementptr i32, i32* %140, i64 4
  %143 = bitcast i32* %142 to <4 x i32>*
  store <4 x i32> %broadcast.splat274, <4 x i32>* %143, align 2
  %index.next266.3 = add i64 %index265, 32
  %144 = getelementptr inbounds i32, i32* %add.ptr.i182, i64 %index.next266.3
  %145 = bitcast i32* %144 to <4 x i32>*
  store <4 x i32> %broadcast.splat274, <4 x i32>* %145, align 2
  %146 = getelementptr i32, i32* %144, i64 4
  %147 = bitcast i32* %146 to <4 x i32>*
  store <4 x i32> %broadcast.splat274, <4 x i32>* %147, align 2
  %index.next266.4 = add i64 %index265, 40
  %148 = getelementptr inbounds i32, i32* %add.ptr.i182, i64 %index.next266.4
  %149 = bitcast i32* %148 to <4 x i32>*
  store <4 x i32> %broadcast.splat274, <4 x i32>* %149, align 2
  %150 = getelementptr i32, i32* %148, i64 4
  %151 = bitcast i32* %150 to <4 x i32>*
  store <4 x i32> %broadcast.splat274, <4 x i32>* %151, align 2
  %index.next266.5 = add i64 %index265, 48
  %152 = getelementptr inbounds i32, i32* %add.ptr.i182, i64 %index.next266.5
  %153 = bitcast i32* %152 to <4 x i32>*
  store <4 x i32> %broadcast.splat274, <4 x i32>* %153, align 2
  %154 = getelementptr i32, i32* %152, i64 4
  %155 = bitcast i32* %154 to <4 x i32>*
  store <4 x i32> %broadcast.splat274, <4 x i32>* %155, align 2
  %index.next266.6 = add i64 %index265, 56
  %156 = getelementptr inbounds i32, i32* %add.ptr.i182, i64 %index.next266.6
  %157 = bitcast i32* %156 to <4 x i32>*
  store <4 x i32> %broadcast.splat274, <4 x i32>* %157, align 2
  %158 = getelementptr i32, i32* %156, i64 4
  %159 = bitcast i32* %158 to <4 x i32>*
  store <4 x i32> %broadcast.splat274, <4 x i32>* %159, align 2
  %index.next266.7 = add i64 %index265, 64
  %160 = icmp eq i64 %index.next266.7, %n.vec262
  br i1 %160, label %middle.block257.unr-lcssa, label %vector.body256, !llvm.loop !30

middle.block257.unr-lcssa:                        ; preds = %vector.body256
  br label %middle.block257

middle.block257:                                  ; preds = %vector.ph264.split, %middle.block257.unr-lcssa
  %cmp.n268 = icmp eq i64 %115, %n.vec262
  br i1 %cmp.n268, label %for.cond5.preheader.i.i, label %for.body.i.i.preheader

for.body.i.i.preheader:                           ; preds = %middle.block257, %min.iters.checked260, %for.body.lr.ph.i.i
  %indvars.iv75.i.i.ph = phi i64 [ 0, %min.iters.checked260 ], [ 0, %for.body.lr.ph.i.i ], [ %n.vec262, %middle.block257 ]
  br label %for.body.i.i

for.body.i.i:                                     ; preds = %for.body.i.i.preheader, %for.body.i.i
  %indvars.iv75.i.i = phi i64 [ %indvars.iv.next76.i.i, %for.body.i.i ], [ %indvars.iv75.i.i.ph, %for.body.i.i.preheader ]
  %arrayidx4.i.i184 = getelementptr inbounds i32, i32* %add.ptr.i182, i64 %indvars.iv75.i.i
  store i32 %DElt.i.sroa.0.0.insert.insert90.i, i32* %arrayidx4.i.i184, align 2
  %indvars.iv.next76.i.i = add nuw nsw i64 %indvars.iv75.i.i, 1
  %lftr.wideiv = trunc i64 %indvars.iv.next76.i.i to i32
  %exitcond = icmp eq i32 %lftr.wideiv, %112
  br i1 %exitcond, label %for.cond5.preheader.i.i.loopexit, label %for.body.i.i, !llvm.loop !31

for.body8.i.i:                                    ; preds = %do.end.i.i, %for.body8.i.preheader.i
  %indvars.iv.i.i = phi i64 [ %indvars.iv.next.i.i, %do.end.i.i ], [ 0, %for.body8.i.preheader.i ]
  %symbol11.i.i = getelementptr inbounds %struct.sortedSymbol_t, %struct.sortedSymbol_t* %add.ptr24.i, i64 %indvars.iv.i.i, i32 0
  %161 = load i8, i8* %symbol11.i.i, align 2, !tbaa !20
  %conv12.i.i = zext i8 %161 to i32
  %weight15.i.i = getelementptr inbounds %struct.sortedSymbol_t, %struct.sortedSymbol_t* %add.ptr24.i, i64 %indvars.iv.i.i, i32 1
  %162 = load i8, i8* %weight15.i.i, align 1, !tbaa !22
  %conv16.i.i = zext i8 %162 to i32
  %sub.i.i = sub i32 %add98, %conv16.i.i
  %sub19.i.i = sub i32 %sub10.i, %sub.i.i
  %shl.i.i = shl i32 1, %sub19.i.i
  %idxprom20.i.i = zext i8 %162 to i64
  %arrayidx21.i.i = getelementptr inbounds [17 x i32], [17 x i32]* %rankVal.i.i, i64 0, i64 %idxprom20.i.i
  %163 = load i32, i32* %arrayidx21.i.i, align 4, !tbaa !2
  %add.i.i = add i32 %shl.i.i, %163
  %shl25.i.i = shl nuw nsw i32 %conv12.i.i, 8
  %add28.i.i = add i32 %sub.i.i, %sub7.i
  %DElt.i.sroa.6.0.insert.ext.i = shl i32 %add28.i.i, 16
  %DElt.i.sroa.6.0.insert.shift.i = and i32 %DElt.i.sroa.6.0.insert.ext.i, 16711680
  %add26.i.i = or i32 %conv.i, %shl25.i.i
  %DElt.i.sroa.6.0.insert.insert.i = or i32 %add26.i.i, %DElt.i.sroa.6.0.insert.shift.i
  %DElt.i.sroa.0.0.insert.insert.i = or i32 %DElt.i.sroa.6.0.insert.insert.i, 33554432
  br label %do.body.i.i

do.body.i.i:                                      ; preds = %do.body.i.i, %for.body8.i.i
  %i22.0.i.i = phi i32 [ %163, %for.body8.i.i ], [ %inc32.i.i, %do.body.i.i ]
  %inc32.i.i = add i32 %i22.0.i.i, 1
  %idxprom33.i.i = zext i32 %i22.0.i.i to i64
  %arrayidx34.i.i183 = getelementptr inbounds i32, i32* %add.ptr.i182, i64 %idxprom33.i.i
  store i32 %DElt.i.sroa.0.0.insert.insert.i, i32* %arrayidx34.i.i183, align 2
  %cmp35.i.i = icmp ult i32 %inc32.i.i, %add.i.i
  br i1 %cmp35.i.i, label %do.body.i.i, label %do.end.i.i

do.end.i.i:                                       ; preds = %do.body.i.i
  store i32 %add.i.i, i32* %arrayidx21.i.i, align 4, !tbaa !2
  %indvars.iv.next.i.i = add nuw nsw i64 %indvars.iv.i.i, 1
  %lftr.wideiv200 = trunc i64 %indvars.iv.next.i.i to i32
  %exitcond201 = icmp eq i32 %lftr.wideiv200, %111
  br i1 %exitcond201, label %HUF_fillDTableX4Level2.exit.i.loopexit, label %for.body8.i.i

HUF_fillDTableX4Level2.exit.i.loopexit:           ; preds = %do.end.i.i
  br label %HUF_fillDTableX4Level2.exit.i

HUF_fillDTableX4Level2.exit.i:                    ; preds = %HUF_fillDTableX4Level2.exit.i.loopexit, %for.cond5.preheader.i.i
  call void @llvm.lifetime.end(i64 68, i8* %105) #11
  %.pre.i = add i32 %shl.i, %108
  br label %if.end36.i

if.else.i:                                        ; preds = %for.body.i
  %add29.i = add i32 %shl.i, %108
  %cmp31101.i = icmp ult i32 %108, %add29.i
  br i1 %cmp31101.i, label %for.body33.lr.ph.i, label %if.end36.i

for.body33.lr.ph.i:                               ; preds = %if.else.i
  %DElt.i.sroa.4.0.insert.ext = shl i32 %sub7.i, 16
  %DElt.i.sroa.4.0.insert.shift = and i32 %DElt.i.sroa.4.0.insert.ext, 16711680
  %DElt.i.sroa.4.0.insert.insert = or i32 %conv.i, %DElt.i.sroa.4.0.insert.shift
  %DElt.i.sroa.0.0.insert.insert = or i32 %DElt.i.sroa.4.0.insert.insert, 16777216
  %164 = zext i32 %108 to i64
  %165 = add i32 %shl.i, -1
  %166 = zext i32 %165 to i64
  %167 = add nuw nsw i64 %166, 1
  %min.iters.check236 = icmp ult i64 %167, 8
  br i1 %min.iters.check236, label %for.body33.i.preheader, label %min.iters.checked237

min.iters.checked237:                             ; preds = %for.body33.lr.ph.i
  %n.vec239 = and i64 %167, 8589934584
  %cmp.zero240 = icmp eq i64 %n.vec239, 0
  %ind.end245 = add nuw nsw i64 %164, %n.vec239
  br i1 %cmp.zero240, label %for.body33.i.preheader, label %vector.ph241

vector.ph241:                                     ; preds = %min.iters.checked237
  %broadcast.splatinsert252 = insertelement <4 x i32> undef, i32 %DElt.i.sroa.0.0.insert.insert, i32 0
  %broadcast.splat253 = shufflevector <4 x i32> %broadcast.splatinsert252, <4 x i32> undef, <4 x i32> zeroinitializer
  %168 = add i32 %shl.i, -1
  %169 = zext i32 %168 to i64
  %170 = add nuw nsw i64 %169, 1
  %171 = and i64 %170, 8589934584
  %172 = add nsw i64 %171, -8
  %173 = lshr exact i64 %172, 3
  %174 = add nuw nsw i64 %173, 1
  %xtraiter277 = and i64 %174, 3
  %lcmp.mod278 = icmp eq i64 %xtraiter277, 0
  br i1 %lcmp.mod278, label %vector.ph241.split, label %vector.body233.prol.preheader

vector.body233.prol.preheader:                    ; preds = %vector.ph241
  br label %vector.body233.prol

vector.body233.prol:                              ; preds = %vector.body233.prol, %vector.body233.prol.preheader
  %index242.prol = phi i64 [ %index.next243.prol, %vector.body233.prol ], [ 0, %vector.body233.prol.preheader ]
  %prol.iter279 = phi i64 [ %prol.iter279.sub, %vector.body233.prol ], [ %xtraiter277, %vector.body233.prol.preheader ]
  %offset.idx247.prol = add i64 %164, %index242.prol
  %175 = getelementptr inbounds i32, i32* %add.ptr2, i64 %offset.idx247.prol
  %176 = bitcast i32* %175 to <4 x i32>*
  store <4 x i32> %broadcast.splat253, <4 x i32>* %176, align 2
  %177 = getelementptr i32, i32* %175, i64 4
  %178 = bitcast i32* %177 to <4 x i32>*
  store <4 x i32> %broadcast.splat253, <4 x i32>* %178, align 2
  %index.next243.prol = add i64 %index242.prol, 8
  %prol.iter279.sub = add i64 %prol.iter279, -1
  %prol.iter279.cmp = icmp eq i64 %prol.iter279.sub, 0
  br i1 %prol.iter279.cmp, label %vector.ph241.split.loopexit, label %vector.body233.prol, !llvm.loop !33

vector.ph241.split.loopexit:                      ; preds = %vector.body233.prol
  %index.next243.prol.lcssa = phi i64 [ %index.next243.prol, %vector.body233.prol ]
  br label %vector.ph241.split

vector.ph241.split:                               ; preds = %vector.ph241, %vector.ph241.split.loopexit
  %index242.unr = phi i64 [ 0, %vector.ph241 ], [ %index.next243.prol.lcssa, %vector.ph241.split.loopexit ]
  %179 = icmp ult i64 %172, 24
  br i1 %179, label %middle.block234, label %vector.ph241.split.split

vector.ph241.split.split:                         ; preds = %vector.ph241.split
  br label %vector.body233

vector.body233:                                   ; preds = %vector.body233, %vector.ph241.split.split
  %index242 = phi i64 [ %index242.unr, %vector.ph241.split.split ], [ %index.next243.3, %vector.body233 ]
  %offset.idx247 = add i64 %164, %index242
  %180 = getelementptr inbounds i32, i32* %add.ptr2, i64 %offset.idx247
  %181 = bitcast i32* %180 to <4 x i32>*
  store <4 x i32> %broadcast.splat253, <4 x i32>* %181, align 2
  %182 = getelementptr i32, i32* %180, i64 4
  %183 = bitcast i32* %182 to <4 x i32>*
  store <4 x i32> %broadcast.splat253, <4 x i32>* %183, align 2
  %index.next243 = add i64 %index242, 8
  %offset.idx247.1 = add i64 %164, %index.next243
  %184 = getelementptr inbounds i32, i32* %add.ptr2, i64 %offset.idx247.1
  %185 = bitcast i32* %184 to <4 x i32>*
  store <4 x i32> %broadcast.splat253, <4 x i32>* %185, align 2
  %186 = getelementptr i32, i32* %184, i64 4
  %187 = bitcast i32* %186 to <4 x i32>*
  store <4 x i32> %broadcast.splat253, <4 x i32>* %187, align 2
  %index.next243.1 = add i64 %index242, 16
  %offset.idx247.2 = add i64 %164, %index.next243.1
  %188 = getelementptr inbounds i32, i32* %add.ptr2, i64 %offset.idx247.2
  %189 = bitcast i32* %188 to <4 x i32>*
  store <4 x i32> %broadcast.splat253, <4 x i32>* %189, align 2
  %190 = getelementptr i32, i32* %188, i64 4
  %191 = bitcast i32* %190 to <4 x i32>*
  store <4 x i32> %broadcast.splat253, <4 x i32>* %191, align 2
  %index.next243.2 = add i64 %index242, 24
  %offset.idx247.3 = add i64 %164, %index.next243.2
  %192 = getelementptr inbounds i32, i32* %add.ptr2, i64 %offset.idx247.3
  %193 = bitcast i32* %192 to <4 x i32>*
  store <4 x i32> %broadcast.splat253, <4 x i32>* %193, align 2
  %194 = getelementptr i32, i32* %192, i64 4
  %195 = bitcast i32* %194 to <4 x i32>*
  store <4 x i32> %broadcast.splat253, <4 x i32>* %195, align 2
  %index.next243.3 = add i64 %index242, 32
  %196 = icmp eq i64 %index.next243.3, %n.vec239
  br i1 %196, label %middle.block234.unr-lcssa, label %vector.body233, !llvm.loop !34

middle.block234.unr-lcssa:                        ; preds = %vector.body233
  br label %middle.block234

middle.block234:                                  ; preds = %vector.ph241.split, %middle.block234.unr-lcssa
  %cmp.n246 = icmp eq i64 %167, %n.vec239
  br i1 %cmp.n246, label %if.end36.i, label %for.body33.i.preheader

for.body33.i.preheader:                           ; preds = %middle.block234, %min.iters.checked237, %for.body33.lr.ph.i
  %indvars.iv.i.ph = phi i64 [ %164, %min.iters.checked237 ], [ %164, %for.body33.lr.ph.i ], [ %ind.end245, %middle.block234 ]
  br label %for.body33.i

for.body33.i:                                     ; preds = %for.body33.i.preheader, %for.body33.i
  %indvars.iv.i = phi i64 [ %indvars.iv.next.i, %for.body33.i ], [ %indvars.iv.i.ph, %for.body33.i.preheader ]
  %arrayidx35.i185 = getelementptr inbounds i32, i32* %add.ptr2, i64 %indvars.iv.i
  store i32 %DElt.i.sroa.0.0.insert.insert, i32* %arrayidx35.i185, align 2
  %indvars.iv.next.i = add nuw nsw i64 %indvars.iv.i, 1
  %lftr.wideiv107.i = trunc i64 %indvars.iv.next.i to i32
  %exitcond108.i = icmp eq i32 %lftr.wideiv107.i, %add29.i
  br i1 %exitcond108.i, label %if.end36.i.loopexit, label %for.body33.i, !llvm.loop !35

if.end36.i.loopexit:                              ; preds = %for.body33.i
  br label %if.end36.i

if.end36.i:                                       ; preds = %if.end36.i.loopexit, %middle.block234, %if.else.i, %HUF_fillDTableX4Level2.exit.i
  %add39.pre-phi.i = phi i32 [ %.pre.i, %HUF_fillDTableX4Level2.exit.i ], [ %add29.i, %if.else.i ], [ %add29.i, %middle.block234 ], [ %add29.i, %if.end36.i.loopexit ]
  store i32 %add39.pre-phi.i, i32* %arrayidx9.i, align 4, !tbaa !2
  %indvars.iv.next110.i = add nuw nsw i64 %indvars.iv109.i, 1
  %lftr.wideiv202 = trunc i64 %indvars.iv.next110.i to i32
  %exitcond203 = icmp eq i32 %lftr.wideiv202, %nextRankStart.0.lcssa
  br i1 %exitcond203, label %HUF_fillDTableX4.exit.loopexit, label %for.body.i

HUF_fillDTableX4.exit.loopexit:                   ; preds = %if.end36.i
  br label %HUF_fillDTableX4.exit

HUF_fillDTableX4.exit:                            ; preds = %HUF_fillDTableX4.exit.loopexit, %for.end94
  call void @llvm.lifetime.end(i64 68, i8* %104) #11
  %tableLog100 = getelementptr inbounds %struct.DTableDesc, %struct.DTableDesc* %tmpcast, i64 0, i32 2
  store i8 %trunc, i8* %tableLog100, align 2, !tbaa !8
  %tableType = getelementptr inbounds %struct.DTableDesc, %struct.DTableDesc* %tmpcast, i64 0, i32 1
  store i8 1, i8* %tableType, align 1, !tbaa !6
  %197 = bitcast i32* %DTable to i8*
  %198 = call i64 @llvm.objectsize.i64.p0i8(i8* %197, i1 false)
  %call101 = call i8* @__memcpy_chk(i8* %197, i8* %7, i64 4, i64 %198) #11
  br label %cleanup

cleanup:                                          ; preds = %if.end9, %if.end, %entry, %HUF_fillDTableX4.exit
  %retval.0 = phi i64 [ %call6, %HUF_fillDTableX4.exit ], [ -16, %entry ], [ %call6, %if.end ], [ -16, %if.end9 ]
  call void @llvm.lifetime.end(i64 4, i8* %7) #11
  call void @llvm.lifetime.end(i64 4, i8* %6) #11
  call void @llvm.lifetime.end(i64 4, i8* %5) #11
  call void @llvm.lifetime.end(i64 1088, i8* %4) #11
  call void @llvm.lifetime.end(i64 72, i8* %3) #11
  call void @llvm.lifetime.end(i64 68, i8* %2) #11
  call void @llvm.lifetime.end(i64 512, i8* %1) #11
  call void @llvm.lifetime.end(i64 256, i8* %0) #11
  ret i64 %retval.0
}

; Function Attrs: nounwind ssp uwtable
define i64 @HUF_decompress1X4_usingDTable(i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize, i32* %DTable) #0 {
entry:
  %bitD.i = alloca %struct.BIT_DStream_t, align 8
  %DTable.val = load i32, i32* %DTable, align 1
  %conv = and i32 %DTable.val, 65280
  %cmp = icmp eq i32 %conv, 256
  br i1 %cmp, label %if.end, label %cleanup

if.end:                                           ; preds = %entry
  %0 = bitcast %struct.BIT_DStream_t* %bitD.i to i8*
  call void @llvm.lifetime.start(i64 32, i8* %0) #11
  %call.i = call fastcc i64 @BIT_initDStream(%struct.BIT_DStream_t* nonnull %bitD.i, i8* %cSrc, i64 %cSrcSize) #11
  %call1.i = call i32 @HUF_isError(i64 %call.i) #11
  %tobool.i = icmp eq i32 %call1.i, 0
  br i1 %tobool.i, label %cleanup.cont.i, label %HUF_decompress1X4_usingDTable_internal.exit

cleanup.cont.i:                                   ; preds = %if.end
  %add.ptr.i = getelementptr inbounds i8, i8* %dst, i64 %dstSize
  %add.ptr2.i = getelementptr inbounds i32, i32* %DTable, i64 1
  %1 = bitcast i32* %add.ptr2.i to %struct.HUF_DEltX4*
  %DTable.val.i = load i32, i32* %DTable, align 1
  %dtd.sroa.3.0.extract.shift.i = lshr i32 %DTable.val.i, 16
  %conv.i = and i32 %dtd.sroa.3.0.extract.shift.i, 255
  call fastcc void @HUF_decodeStreamX4(i8* %dst, %struct.BIT_DStream_t* nonnull %bitD.i, i8* %add.ptr.i, %struct.HUF_DEltX4* %1, i32 %conv.i) #11
  %bitD.idx.i = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD.i, i64 0, i32 1
  %bitD.idx.val.i = load i32, i32* %bitD.idx.i, align 8, !tbaa !13
  %bitD.idx21.i = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD.i, i64 0, i32 2
  %bitD.idx21.val.i = load i8*, i8** %bitD.idx21.i, align 8, !tbaa !17
  %bitD.idx22.i = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD.i, i64 0, i32 3
  %bitD.idx22.val.i = load i8*, i8** %bitD.idx22.i, align 8, !tbaa !18
  %cmp.i.i = icmp eq i8* %bitD.idx21.val.i, %bitD.idx22.val.i
  %cmp1.i.i = icmp eq i32 %bitD.idx.val.i, 64
  %tobool6.demorgan.i = and i1 %cmp1.i.i, %cmp.i.i
  %.dstSize.i = select i1 %tobool6.demorgan.i, i64 %dstSize, i64 -14
  br label %HUF_decompress1X4_usingDTable_internal.exit

HUF_decompress1X4_usingDTable_internal.exit:      ; preds = %if.end, %cleanup.cont.i
  %retval.1.i = phi i64 [ %call.i, %if.end ], [ %.dstSize.i, %cleanup.cont.i ]
  call void @llvm.lifetime.end(i64 32, i8* %0) #11
  br label %cleanup

cleanup:                                          ; preds = %entry, %HUF_decompress1X4_usingDTable_internal.exit
  %retval.0 = phi i64 [ %retval.1.i, %HUF_decompress1X4_usingDTable_internal.exit ], [ -1, %entry ]
  ret i64 %retval.0
}

; Function Attrs: nounwind ssp uwtable
define i64 @HUF_decompress1X4_DCtx(i32* %DCtx, i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize) #0 {
entry:
  %bitD.i = alloca %struct.BIT_DStream_t, align 8
  %call = tail call i64 @HUF_readDTableX4(i32* %DCtx, i8* %cSrc, i64 %cSrcSize)
  %call1 = tail call i32 @HUF_isError(i64 %call) #11
  %tobool = icmp eq i32 %call1, 0
  br i1 %tobool, label %if.end, label %cleanup

if.end:                                           ; preds = %entry
  %cmp = icmp ult i64 %call, %cSrcSize
  br i1 %cmp, label %if.end3, label %cleanup

if.end3:                                          ; preds = %if.end
  %add.ptr = getelementptr inbounds i8, i8* %cSrc, i64 %call
  %sub = sub i64 %cSrcSize, %call
  %0 = bitcast %struct.BIT_DStream_t* %bitD.i to i8*
  call void @llvm.lifetime.start(i64 32, i8* %0) #11
  %call.i = call fastcc i64 @BIT_initDStream(%struct.BIT_DStream_t* nonnull %bitD.i, i8* %add.ptr, i64 %sub) #11
  %call1.i = call i32 @HUF_isError(i64 %call.i) #11
  %tobool.i = icmp eq i32 %call1.i, 0
  br i1 %tobool.i, label %cleanup.cont.i, label %HUF_decompress1X4_usingDTable_internal.exit

cleanup.cont.i:                                   ; preds = %if.end3
  %add.ptr.i = getelementptr inbounds i8, i8* %dst, i64 %dstSize
  %add.ptr2.i = getelementptr inbounds i32, i32* %DCtx, i64 1
  %1 = bitcast i32* %add.ptr2.i to %struct.HUF_DEltX4*
  %DTable.val.i = load i32, i32* %DCtx, align 1
  %dtd.sroa.3.0.extract.shift.i = lshr i32 %DTable.val.i, 16
  %conv.i = and i32 %dtd.sroa.3.0.extract.shift.i, 255
  call fastcc void @HUF_decodeStreamX4(i8* %dst, %struct.BIT_DStream_t* nonnull %bitD.i, i8* %add.ptr.i, %struct.HUF_DEltX4* %1, i32 %conv.i) #11
  %bitD.idx.i = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD.i, i64 0, i32 1
  %bitD.idx.val.i = load i32, i32* %bitD.idx.i, align 8, !tbaa !13
  %bitD.idx21.i = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD.i, i64 0, i32 2
  %bitD.idx21.val.i = load i8*, i8** %bitD.idx21.i, align 8, !tbaa !17
  %bitD.idx22.i = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD.i, i64 0, i32 3
  %bitD.idx22.val.i = load i8*, i8** %bitD.idx22.i, align 8, !tbaa !18
  %cmp.i.i = icmp eq i8* %bitD.idx21.val.i, %bitD.idx22.val.i
  %cmp1.i.i = icmp eq i32 %bitD.idx.val.i, 64
  %tobool6.demorgan.i = and i1 %cmp1.i.i, %cmp.i.i
  %.dstSize.i = select i1 %tobool6.demorgan.i, i64 %dstSize, i64 -14
  br label %HUF_decompress1X4_usingDTable_internal.exit

HUF_decompress1X4_usingDTable_internal.exit:      ; preds = %if.end3, %cleanup.cont.i
  %retval.1.i = phi i64 [ %call.i, %if.end3 ], [ %.dstSize.i, %cleanup.cont.i ]
  call void @llvm.lifetime.end(i64 32, i8* %0) #11
  br label %cleanup

cleanup:                                          ; preds = %if.end, %entry, %HUF_decompress1X4_usingDTable_internal.exit
  %retval.0 = phi i64 [ %retval.1.i, %HUF_decompress1X4_usingDTable_internal.exit ], [ %call, %entry ], [ -13, %if.end ]
  ret i64 %retval.0
}

; Function Attrs: nounwind ssp uwtable
define i64 @HUF_decompress1X4(i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize) #0 {
entry:
  %DTable = alloca [4097 x i32], align 16
  %0 = bitcast [4097 x i32]* %DTable to i8*
  call void @llvm.lifetime.start(i64 16388, i8* %0) #11
  call void @llvm.memset.p0i8.i64(i8* %0, i8 0, i64 16388, i32 16, i1 false)
  %1 = getelementptr [4097 x i32], [4097 x i32]* %DTable, i64 0, i64 0
  store i32 201326604, i32* %1, align 16
  %call = call i64 @HUF_decompress1X4_DCtx(i32* %1, i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize)
  call void @llvm.lifetime.end(i64 16388, i8* %0) #11
  ret i64 %call
}

; Function Attrs: nounwind ssp uwtable
define i64 @HUF_decompress4X4_usingDTable(i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize, i32* %DTable) #0 {
entry:
  %DTable.val = load i32, i32* %DTable, align 1
  %conv = and i32 %DTable.val, 65280
  %cmp = icmp eq i32 %conv, 256
  br i1 %cmp, label %if.end, label %cleanup

if.end:                                           ; preds = %entry
  %call2 = tail call fastcc i64 @HUF_decompress4X4_usingDTable_internal(i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize, i32* nonnull %DTable)
  br label %cleanup

cleanup:                                          ; preds = %entry, %if.end
  %retval.0 = phi i64 [ %call2, %if.end ], [ -1, %entry ]
  ret i64 %retval.0
}

; Function Attrs: nounwind ssp uwtable
define internal fastcc i64 @HUF_decompress4X4_usingDTable_internal(i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize, i32* %DTable) unnamed_addr #0 {
entry:
  %bitD1 = alloca %struct.BIT_DStream_t, align 8
  %bitD2 = alloca %struct.BIT_DStream_t, align 8
  %bitD3 = alloca %struct.BIT_DStream_t, align 8
  %bitD4 = alloca %struct.BIT_DStream_t, align 8
  %cmp = icmp ult i64 %cSrcSize, 10
  br i1 %cmp, label %return, label %if.end

if.end:                                           ; preds = %entry
  %add.ptr = getelementptr inbounds i8, i8* %dst, i64 %dstSize
  %add.ptr1 = getelementptr inbounds i32, i32* %DTable, i64 1
  %0 = bitcast i32* %add.ptr1 to %struct.HUF_DEltX4*
  %1 = bitcast %struct.BIT_DStream_t* %bitD1 to i8*
  call void @llvm.lifetime.start(i64 32, i8* %1) #11
  %2 = bitcast %struct.BIT_DStream_t* %bitD2 to i8*
  call void @llvm.lifetime.start(i64 32, i8* %2) #11
  %3 = bitcast %struct.BIT_DStream_t* %bitD3 to i8*
  call void @llvm.lifetime.start(i64 32, i8* %3) #11
  %4 = bitcast %struct.BIT_DStream_t* %bitD4 to i8*
  call void @llvm.lifetime.start(i64 32, i8* %4) #11
  %val.0..sroa_cast1.i.i = bitcast i8* %cSrc to i16*
  %val.0.copyload.i.i = load i16, i16* %val.0..sroa_cast1.i.i, align 1
  %conv = zext i16 %val.0.copyload.i.i to i64
  %add.ptr2 = getelementptr inbounds i8, i8* %cSrc, i64 2
  %val.0..sroa_cast1.i.i402 = bitcast i8* %add.ptr2 to i16*
  %val.0.copyload.i.i403 = load i16, i16* %val.0..sroa_cast1.i.i402, align 1
  %conv4 = zext i16 %val.0.copyload.i.i403 to i64
  %add.ptr5 = getelementptr inbounds i8, i8* %cSrc, i64 4
  %val.0..sroa_cast1.i.i404 = bitcast i8* %add.ptr5 to i16*
  %val.0.copyload.i.i405 = load i16, i16* %val.0..sroa_cast1.i.i404, align 1
  %conv7 = zext i16 %val.0.copyload.i.i405 to i64
  %add.neg = add i64 %cSrcSize, -6
  %add8.neg = sub i64 %add.neg, %conv
  %add9.neg = sub i64 %add8.neg, %conv4
  %sub = sub i64 %add9.neg, %conv7
  %add.ptr10 = getelementptr inbounds i8, i8* %cSrc, i64 6
  %add.ptr11 = getelementptr inbounds i8, i8* %add.ptr10, i64 %conv
  %add.ptr12 = getelementptr inbounds i8, i8* %add.ptr11, i64 %conv4
  %add.ptr13 = getelementptr inbounds i8, i8* %add.ptr12, i64 %conv7
  %add14 = add i64 %dstSize, 3
  %div = lshr i64 %add14, 2
  %add.ptr15 = getelementptr inbounds i8, i8* %dst, i64 %div
  %add.ptr16 = getelementptr inbounds i8, i8* %add.ptr15, i64 %div
  %add.ptr17 = getelementptr inbounds i8, i8* %add.ptr16, i64 %div
  %DTable.val = load i32, i32* %DTable, align 1
  %dtd.sroa.3.0.extract.shift = lshr i32 %DTable.val, 16
  %conv19 = and i32 %dtd.sroa.3.0.extract.shift, 255
  %cmp20 = icmp ugt i64 %sub, %cSrcSize
  br i1 %cmp20, label %cleanup201, label %if.end23

if.end23:                                         ; preds = %if.end
  %call24 = call fastcc i64 @BIT_initDStream(%struct.BIT_DStream_t* nonnull %bitD1, i8* %add.ptr10, i64 %conv)
  %call25 = call i32 @HUF_isError(i64 %call24) #11
  %tobool = icmp eq i32 %call25, 0
  br i1 %tobool, label %cleanup.cont, label %cleanup201

cleanup.cont:                                     ; preds = %if.end23
  %call29 = call fastcc i64 @BIT_initDStream(%struct.BIT_DStream_t* nonnull %bitD2, i8* %add.ptr11, i64 %conv4)
  %call30 = call i32 @HUF_isError(i64 %call29) #11
  %tobool31 = icmp eq i32 %call30, 0
  br i1 %tobool31, label %cleanup.cont36, label %cleanup201

cleanup.cont36:                                   ; preds = %cleanup.cont
  %call38 = call fastcc i64 @BIT_initDStream(%struct.BIT_DStream_t* nonnull %bitD3, i8* %add.ptr12, i64 %conv7)
  %call39 = call i32 @HUF_isError(i64 %call38) #11
  %tobool40 = icmp eq i32 %call39, 0
  br i1 %tobool40, label %cleanup.cont45, label %cleanup201

cleanup.cont45:                                   ; preds = %cleanup.cont36
  %call47 = call fastcc i64 @BIT_initDStream(%struct.BIT_DStream_t* nonnull %bitD4, i8* %add.ptr13, i64 %sub)
  %call48 = call i32 @HUF_isError(i64 %call47) #11
  %tobool49 = icmp eq i32 %call48, 0
  br i1 %tobool49, label %cleanup.cont54, label %cleanup201

cleanup.cont54:                                   ; preds = %cleanup.cont45
  %call55 = call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* nonnull %bitD1)
  %call56 = call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* nonnull %bitD2)
  %or = or i32 %call56, %call55
  %call57 = call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* nonnull %bitD3)
  %or58 = or i32 %or, %call57
  %call59 = call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* nonnull %bitD4)
  %or60 = or i32 %or58, %call59
  %cmp61409 = icmp eq i32 %or60, 0
  %add.ptr63 = getelementptr inbounds i8, i8* %add.ptr, i64 -7
  %cmp64410 = icmp ult i8* %add.ptr17, %add.ptr63
  %and376411 = and i1 %cmp61409, %cmp64410
  br i1 %and376411, label %for.body.preheader, label %for.end

for.body.preheader:                               ; preds = %cleanup.cont54
  br label %for.body

for.body:                                         ; preds = %for.body.preheader, %for.body
  %op1.0415 = phi i8* [ %add.ptr155, %for.body ], [ %dst, %for.body.preheader ]
  %op2.0414 = phi i8* [ %add.ptr158, %for.body ], [ %add.ptr15, %for.body.preheader ]
  %op3.0413 = phi i8* [ %add.ptr161, %for.body ], [ %add.ptr16, %for.body.preheader ]
  %op4.0412 = phi i8* [ %add.ptr164, %for.body ], [ %add.ptr17, %for.body.preheader ]
  %call70 = call fastcc i32 @HUF_decodeSymbolX4(i8* %op1.0415, %struct.BIT_DStream_t* nonnull %bitD1, %struct.HUF_DEltX4* %0, i32 %conv19)
  %idx.ext = zext i32 %call70 to i64
  %add.ptr71 = getelementptr inbounds i8, i8* %op1.0415, i64 %idx.ext
  %call76 = call fastcc i32 @HUF_decodeSymbolX4(i8* %op2.0414, %struct.BIT_DStream_t* nonnull %bitD2, %struct.HUF_DEltX4* %0, i32 %conv19)
  %idx.ext77 = zext i32 %call76 to i64
  %add.ptr78 = getelementptr inbounds i8, i8* %op2.0414, i64 %idx.ext77
  %call83 = call fastcc i32 @HUF_decodeSymbolX4(i8* %op3.0413, %struct.BIT_DStream_t* nonnull %bitD3, %struct.HUF_DEltX4* %0, i32 %conv19)
  %idx.ext84 = zext i32 %call83 to i64
  %add.ptr85 = getelementptr inbounds i8, i8* %op3.0413, i64 %idx.ext84
  %call90 = call fastcc i32 @HUF_decodeSymbolX4(i8* %op4.0412, %struct.BIT_DStream_t* nonnull %bitD4, %struct.HUF_DEltX4* %0, i32 %conv19)
  %idx.ext91 = zext i32 %call90 to i64
  %add.ptr92 = getelementptr inbounds i8, i8* %op4.0412, i64 %idx.ext91
  %call97 = call fastcc i32 @HUF_decodeSymbolX4(i8* %add.ptr71, %struct.BIT_DStream_t* nonnull %bitD1, %struct.HUF_DEltX4* %0, i32 %conv19)
  %idx.ext98 = zext i32 %call97 to i64
  %add.ptr99 = getelementptr inbounds i8, i8* %add.ptr71, i64 %idx.ext98
  %call105 = call fastcc i32 @HUF_decodeSymbolX4(i8* %add.ptr78, %struct.BIT_DStream_t* nonnull %bitD2, %struct.HUF_DEltX4* %0, i32 %conv19)
  %idx.ext106 = zext i32 %call105 to i64
  %add.ptr107 = getelementptr inbounds i8, i8* %add.ptr78, i64 %idx.ext106
  %call113 = call fastcc i32 @HUF_decodeSymbolX4(i8* %add.ptr85, %struct.BIT_DStream_t* nonnull %bitD3, %struct.HUF_DEltX4* %0, i32 %conv19)
  %idx.ext114 = zext i32 %call113 to i64
  %add.ptr115 = getelementptr inbounds i8, i8* %add.ptr85, i64 %idx.ext114
  %call121 = call fastcc i32 @HUF_decodeSymbolX4(i8* %add.ptr92, %struct.BIT_DStream_t* nonnull %bitD4, %struct.HUF_DEltX4* %0, i32 %conv19)
  %idx.ext122 = zext i32 %call121 to i64
  %add.ptr123 = getelementptr inbounds i8, i8* %add.ptr92, i64 %idx.ext122
  %call128 = call fastcc i32 @HUF_decodeSymbolX4(i8* %add.ptr99, %struct.BIT_DStream_t* nonnull %bitD1, %struct.HUF_DEltX4* %0, i32 %conv19)
  %idx.ext129 = zext i32 %call128 to i64
  %add.ptr130 = getelementptr inbounds i8, i8* %add.ptr99, i64 %idx.ext129
  %call135 = call fastcc i32 @HUF_decodeSymbolX4(i8* %add.ptr107, %struct.BIT_DStream_t* nonnull %bitD2, %struct.HUF_DEltX4* %0, i32 %conv19)
  %idx.ext136 = zext i32 %call135 to i64
  %add.ptr137 = getelementptr inbounds i8, i8* %add.ptr107, i64 %idx.ext136
  %call142 = call fastcc i32 @HUF_decodeSymbolX4(i8* %add.ptr115, %struct.BIT_DStream_t* nonnull %bitD3, %struct.HUF_DEltX4* %0, i32 %conv19)
  %idx.ext143 = zext i32 %call142 to i64
  %add.ptr144 = getelementptr inbounds i8, i8* %add.ptr115, i64 %idx.ext143
  %call149 = call fastcc i32 @HUF_decodeSymbolX4(i8* %add.ptr123, %struct.BIT_DStream_t* nonnull %bitD4, %struct.HUF_DEltX4* %0, i32 %conv19)
  %idx.ext150 = zext i32 %call149 to i64
  %add.ptr151 = getelementptr inbounds i8, i8* %add.ptr123, i64 %idx.ext150
  %call153 = call fastcc i32 @HUF_decodeSymbolX4(i8* %add.ptr130, %struct.BIT_DStream_t* nonnull %bitD1, %struct.HUF_DEltX4* %0, i32 %conv19)
  %idx.ext154 = zext i32 %call153 to i64
  %add.ptr155 = getelementptr inbounds i8, i8* %add.ptr130, i64 %idx.ext154
  %call156 = call fastcc i32 @HUF_decodeSymbolX4(i8* %add.ptr137, %struct.BIT_DStream_t* nonnull %bitD2, %struct.HUF_DEltX4* %0, i32 %conv19)
  %idx.ext157 = zext i32 %call156 to i64
  %add.ptr158 = getelementptr inbounds i8, i8* %add.ptr137, i64 %idx.ext157
  %call159 = call fastcc i32 @HUF_decodeSymbolX4(i8* %add.ptr144, %struct.BIT_DStream_t* nonnull %bitD3, %struct.HUF_DEltX4* %0, i32 %conv19)
  %idx.ext160 = zext i32 %call159 to i64
  %add.ptr161 = getelementptr inbounds i8, i8* %add.ptr144, i64 %idx.ext160
  %call162 = call fastcc i32 @HUF_decodeSymbolX4(i8* %add.ptr151, %struct.BIT_DStream_t* nonnull %bitD4, %struct.HUF_DEltX4* %0, i32 %conv19)
  %idx.ext163 = zext i32 %call162 to i64
  %add.ptr164 = getelementptr inbounds i8, i8* %add.ptr151, i64 %idx.ext163
  %call165 = call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* nonnull %bitD1)
  %call166 = call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* nonnull %bitD2)
  %or167 = or i32 %call166, %call165
  %call168 = call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* nonnull %bitD3)
  %or169 = or i32 %or167, %call168
  %call170 = call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* nonnull %bitD4)
  %or171 = or i32 %or169, %call170
  %cmp61 = icmp eq i32 %or171, 0
  %cmp64 = icmp ult i8* %add.ptr164, %add.ptr63
  %and376 = and i1 %cmp61, %cmp64
  br i1 %and376, label %for.body, label %for.end.loopexit

for.end.loopexit:                                 ; preds = %for.body
  %add.ptr164.lcssa = phi i8* [ %add.ptr164, %for.body ]
  %add.ptr161.lcssa = phi i8* [ %add.ptr161, %for.body ]
  %add.ptr158.lcssa = phi i8* [ %add.ptr158, %for.body ]
  %add.ptr155.lcssa = phi i8* [ %add.ptr155, %for.body ]
  br label %for.end

for.end:                                          ; preds = %for.end.loopexit, %cleanup.cont54
  %op1.0.lcssa = phi i8* [ %dst, %cleanup.cont54 ], [ %add.ptr155.lcssa, %for.end.loopexit ]
  %op2.0.lcssa = phi i8* [ %add.ptr15, %cleanup.cont54 ], [ %add.ptr158.lcssa, %for.end.loopexit ]
  %op3.0.lcssa = phi i8* [ %add.ptr16, %cleanup.cont54 ], [ %add.ptr161.lcssa, %for.end.loopexit ]
  %op4.0.lcssa = phi i8* [ %add.ptr17, %cleanup.cont54 ], [ %add.ptr164.lcssa, %for.end.loopexit ]
  %cmp172 = icmp ugt i8* %op1.0.lcssa, %add.ptr15
  %cmp176 = icmp ugt i8* %op2.0.lcssa, %add.ptr16
  %or.cond = or i1 %cmp176, %cmp172
  %cmp180 = icmp ugt i8* %op3.0.lcssa, %add.ptr17
  %or.cond381 = or i1 %cmp180, %or.cond
  br i1 %or.cond381, label %cleanup201, label %if.end183

if.end183:                                        ; preds = %for.end
  call fastcc void @HUF_decodeStreamX4(i8* %op1.0.lcssa, %struct.BIT_DStream_t* nonnull %bitD1, i8* %add.ptr15, %struct.HUF_DEltX4* %0, i32 %conv19)
  call fastcc void @HUF_decodeStreamX4(i8* %op2.0.lcssa, %struct.BIT_DStream_t* nonnull %bitD2, i8* %add.ptr16, %struct.HUF_DEltX4* %0, i32 %conv19)
  call fastcc void @HUF_decodeStreamX4(i8* %op3.0.lcssa, %struct.BIT_DStream_t* nonnull %bitD3, i8* %add.ptr17, %struct.HUF_DEltX4* %0, i32 %conv19)
  call fastcc void @HUF_decodeStreamX4(i8* %op4.0.lcssa, %struct.BIT_DStream_t* nonnull %bitD4, i8* %add.ptr, %struct.HUF_DEltX4* %0, i32 %conv19)
  %bitD1.idx = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD1, i64 0, i32 1
  %bitD1.idx.val = load i32, i32* %bitD1.idx, align 8, !tbaa !13
  %bitD1.idx388 = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD1, i64 0, i32 2
  %bitD1.idx388.val = load i8*, i8** %bitD1.idx388, align 8, !tbaa !17
  %bitD1.idx389 = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD1, i64 0, i32 3
  %bitD1.idx389.val = load i8*, i8** %bitD1.idx389, align 8, !tbaa !18
  %cmp.i398 = icmp eq i8* %bitD1.idx388.val, %bitD1.idx389.val
  %cmp1.i399 = icmp eq i32 %bitD1.idx.val, 64
  %cmp1..i400 = and i1 %cmp1.i399, %cmp.i398
  %bitD2.idx = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD2, i64 0, i32 1
  %bitD2.idx.val = load i32, i32* %bitD2.idx, align 8, !tbaa !13
  %bitD2.idx386 = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD2, i64 0, i32 2
  %bitD2.idx386.val = load i8*, i8** %bitD2.idx386, align 8, !tbaa !17
  %bitD2.idx387 = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD2, i64 0, i32 3
  %bitD2.idx387.val = load i8*, i8** %bitD2.idx387, align 8, !tbaa !18
  %cmp.i394 = icmp eq i8* %bitD2.idx386.val, %bitD2.idx387.val
  %cmp1.i395 = icmp eq i32 %bitD2.idx.val, 64
  %cmp1..i396 = and i1 %cmp1.i395, %cmp.i394
  %and190406 = and i1 %cmp1..i400, %cmp1..i396
  %bitD3.idx = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD3, i64 0, i32 1
  %bitD3.idx.val = load i32, i32* %bitD3.idx, align 8, !tbaa !13
  %bitD3.idx384 = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD3, i64 0, i32 2
  %bitD3.idx384.val = load i8*, i8** %bitD3.idx384, align 8, !tbaa !17
  %bitD3.idx385 = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD3, i64 0, i32 3
  %bitD3.idx385.val = load i8*, i8** %bitD3.idx385, align 8, !tbaa !18
  %cmp.i390 = icmp eq i8* %bitD3.idx384.val, %bitD3.idx385.val
  %cmp1.i391 = icmp eq i32 %bitD3.idx.val, 64
  %cmp1..i392 = and i1 %cmp1.i391, %cmp.i390
  %and192407 = and i1 %and190406, %cmp1..i392
  %bitD4.idx = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD4, i64 0, i32 1
  %bitD4.idx.val = load i32, i32* %bitD4.idx, align 8, !tbaa !13
  %bitD4.idx382 = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD4, i64 0, i32 2
  %bitD4.idx382.val = load i8*, i8** %bitD4.idx382, align 8, !tbaa !17
  %bitD4.idx383 = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD4, i64 0, i32 3
  %bitD4.idx383.val = load i8*, i8** %bitD4.idx383, align 8, !tbaa !18
  %cmp.i = icmp eq i8* %bitD4.idx382.val, %bitD4.idx383.val
  %cmp1.i = icmp eq i32 %bitD4.idx.val, 64
  %cmp1..i = and i1 %cmp1.i, %cmp.i
  %and194408 = and i1 %and192407, %cmp1..i
  %dstSize..call24.call29.call38.call47 = select i1 %and194408, i64 %dstSize, i64 -14
  br label %cleanup201

cleanup201:                                       ; preds = %if.end183, %for.end, %if.end, %cleanup.cont45, %cleanup.cont36, %cleanup.cont, %if.end23
  %retval.5 = phi i64 [ %call47, %cleanup.cont45 ], [ %call38, %cleanup.cont36 ], [ %call29, %cleanup.cont ], [ %call24, %if.end23 ], [ -14, %if.end ], [ -14, %for.end ], [ %dstSize..call24.call29.call38.call47, %if.end183 ]
  call void @llvm.lifetime.end(i64 32, i8* %4) #11
  call void @llvm.lifetime.end(i64 32, i8* %3) #11
  call void @llvm.lifetime.end(i64 32, i8* %2) #11
  call void @llvm.lifetime.end(i64 32, i8* %1) #11
  br label %return

return:                                           ; preds = %entry, %cleanup201
  %retval.6 = phi i64 [ %retval.5, %cleanup201 ], [ -14, %entry ]
  ret i64 %retval.6
}

; Function Attrs: nounwind ssp uwtable
define i64 @HUF_decompress4X4_DCtx(i32* %dctx, i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize) #0 {
entry:
  %call = tail call i64 @HUF_readDTableX4(i32* %dctx, i8* %cSrc, i64 %cSrcSize)
  %call1 = tail call i32 @HUF_isError(i64 %call) #11
  %tobool = icmp eq i32 %call1, 0
  br i1 %tobool, label %if.end, label %cleanup

if.end:                                           ; preds = %entry
  %cmp = icmp ult i64 %call, %cSrcSize
  br i1 %cmp, label %if.end3, label %cleanup

if.end3:                                          ; preds = %if.end
  %add.ptr = getelementptr inbounds i8, i8* %cSrc, i64 %call
  %sub = sub i64 %cSrcSize, %call
  %call4 = tail call fastcc i64 @HUF_decompress4X4_usingDTable_internal(i8* %dst, i64 %dstSize, i8* %add.ptr, i64 %sub, i32* %dctx)
  br label %cleanup

cleanup:                                          ; preds = %if.end, %entry, %if.end3
  %retval.0 = phi i64 [ %call4, %if.end3 ], [ %call, %entry ], [ -13, %if.end ]
  ret i64 %retval.0
}

; Function Attrs: nounwind ssp uwtable
define i64 @HUF_decompress4X4(i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize) #0 {
entry:
  %DTable = alloca [4097 x i32], align 16
  %0 = bitcast [4097 x i32]* %DTable to i8*
  call void @llvm.lifetime.start(i64 16388, i8* %0) #11
  call void @llvm.memset.p0i8.i64(i8* %0, i8 0, i64 16388, i32 16, i1 false)
  %1 = getelementptr [4097 x i32], [4097 x i32]* %DTable, i64 0, i64 0
  store i32 201326604, i32* %1, align 16
  %call.i = call i64 @HUF_readDTableX4(i32* %1, i8* %cSrc, i64 %cSrcSize) #11
  %call1.i = call i32 @HUF_isError(i64 %call.i) #11
  %tobool.i = icmp eq i32 %call1.i, 0
  br i1 %tobool.i, label %if.end.i, label %HUF_decompress4X4_DCtx.exit

if.end.i:                                         ; preds = %entry
  %cmp.i = icmp ult i64 %call.i, %cSrcSize
  br i1 %cmp.i, label %if.end3.i, label %HUF_decompress4X4_DCtx.exit

if.end3.i:                                        ; preds = %if.end.i
  %add.ptr.i = getelementptr inbounds i8, i8* %cSrc, i64 %call.i
  %sub.i = sub i64 %cSrcSize, %call.i
  %call4.i = call fastcc i64 @HUF_decompress4X4_usingDTable_internal(i8* %dst, i64 %dstSize, i8* %add.ptr.i, i64 %sub.i, i32* %1) #11
  br label %HUF_decompress4X4_DCtx.exit

HUF_decompress4X4_DCtx.exit:                      ; preds = %entry, %if.end.i, %if.end3.i
  %retval.0.i = phi i64 [ %call4.i, %if.end3.i ], [ %call.i, %entry ], [ -13, %if.end.i ]
  call void @llvm.lifetime.end(i64 16388, i8* %0) #11
  ret i64 %retval.0.i
}

; Function Attrs: nounwind ssp uwtable
define i64 @HUF_decompress1X_usingDTable(i8* %dst, i64 %maxDstSize, i8* %cSrc, i64 %cSrcSize, i32* %DTable) #0 {
entry:
  %bitD.i6 = alloca %struct.BIT_DStream_t, align 8
  %bitD.i = alloca %struct.BIT_DStream_t, align 8
  %DTable.val = load i32, i32* %DTable, align 1
  %conv = and i32 %DTable.val, 65280
  %tobool = icmp eq i32 %conv, 0
  br i1 %tobool, label %cond.false, label %cond.true

cond.true:                                        ; preds = %entry
  %0 = bitcast %struct.BIT_DStream_t* %bitD.i6 to i8*
  call void @llvm.lifetime.start(i64 32, i8* %0) #11
  %call.i = call fastcc i64 @BIT_initDStream(%struct.BIT_DStream_t* nonnull %bitD.i6, i8* %cSrc, i64 %cSrcSize) #11
  %call1.i = call i32 @HUF_isError(i64 %call.i) #11
  %tobool.i7 = icmp eq i32 %call1.i, 0
  br i1 %tobool.i7, label %cleanup.cont.i18, label %HUF_decompress1X4_usingDTable_internal.exit

cleanup.cont.i18:                                 ; preds = %cond.true
  %add.ptr.i8 = getelementptr inbounds i8, i8* %dst, i64 %maxDstSize
  %add.ptr2.i = getelementptr inbounds i32, i32* %DTable, i64 1
  %1 = bitcast i32* %add.ptr2.i to %struct.HUF_DEltX4*
  %DTable.val.i9 = load i32, i32* %DTable, align 1
  %dtd.sroa.3.0.extract.shift.i10 = lshr i32 %DTable.val.i9, 16
  %conv.i11 = and i32 %dtd.sroa.3.0.extract.shift.i10, 255
  call fastcc void @HUF_decodeStreamX4(i8* %dst, %struct.BIT_DStream_t* nonnull %bitD.i6, i8* %add.ptr.i8, %struct.HUF_DEltX4* %1, i32 %conv.i11) #11
  %bitD.idx.i12 = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD.i6, i64 0, i32 1
  %bitD.idx.val.i13 = load i32, i32* %bitD.idx.i12, align 8, !tbaa !13
  %bitD.idx21.i = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD.i6, i64 0, i32 2
  %bitD.idx21.val.i = load i8*, i8** %bitD.idx21.i, align 8, !tbaa !17
  %bitD.idx22.i = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD.i6, i64 0, i32 3
  %bitD.idx22.val.i = load i8*, i8** %bitD.idx22.i, align 8, !tbaa !18
  %cmp.i.i14 = icmp eq i8* %bitD.idx21.val.i, %bitD.idx22.val.i
  %cmp1.i.i15 = icmp eq i32 %bitD.idx.val.i13, 64
  %tobool6.demorgan.i16 = and i1 %cmp1.i.i15, %cmp.i.i14
  %.dstSize.i17 = select i1 %tobool6.demorgan.i16, i64 %maxDstSize, i64 -14
  br label %HUF_decompress1X4_usingDTable_internal.exit

HUF_decompress1X4_usingDTable_internal.exit:      ; preds = %cond.true, %cleanup.cont.i18
  %retval.1.i19 = phi i64 [ %call.i, %cond.true ], [ %.dstSize.i17, %cleanup.cont.i18 ]
  call void @llvm.lifetime.end(i64 32, i8* %0) #11
  br label %cond.end

cond.false:                                       ; preds = %entry
  %2 = bitcast %struct.BIT_DStream_t* %bitD.i to i8*
  call void @llvm.lifetime.start(i64 32, i8* %2) #11
  %call2.i = call fastcc i64 @BIT_initDStream(%struct.BIT_DStream_t* nonnull %bitD.i, i8* %cSrc, i64 %cSrcSize) #11
  %call3.i = call i32 @HUF_isError(i64 %call2.i) #11
  %tobool.i = icmp eq i32 %call3.i, 0
  br i1 %tobool.i, label %cleanup.cont.i, label %HUF_decompress1X2_usingDTable_internal.exit

cleanup.cont.i:                                   ; preds = %cond.false
  %dtd.sroa.3.0.extract.shift.i = lshr i32 %DTable.val, 16
  %conv.i = and i32 %dtd.sroa.3.0.extract.shift.i, 255
  %add.ptr1.i = getelementptr inbounds i32, i32* %DTable, i64 1
  %3 = bitcast i32* %add.ptr1.i to %struct.HUF_DEltX2*
  %add.ptr.i = getelementptr inbounds i8, i8* %dst, i64 %maxDstSize
  call fastcc void @HUF_decodeStreamX2(i8* %dst, %struct.BIT_DStream_t* nonnull %bitD.i, i8* %add.ptr.i, %struct.HUF_DEltX2* %3, i32 %conv.i) #11
  %bitD.idx.i = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD.i, i64 0, i32 1
  %bitD.idx.val.i = load i32, i32* %bitD.idx.i, align 8, !tbaa !13
  %bitD.idx28.i = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD.i, i64 0, i32 2
  %bitD.idx28.val.i = load i8*, i8** %bitD.idx28.i, align 8, !tbaa !17
  %bitD.idx29.i = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD.i, i64 0, i32 3
  %bitD.idx29.val.i = load i8*, i8** %bitD.idx29.i, align 8, !tbaa !18
  %cmp.i.i = icmp eq i8* %bitD.idx28.val.i, %bitD.idx29.val.i
  %cmp1.i.i = icmp eq i32 %bitD.idx.val.i, 64
  %tobool6.demorgan.i = and i1 %cmp1.i.i, %cmp.i.i
  %.dstSize.i = select i1 %tobool6.demorgan.i, i64 %maxDstSize, i64 -14
  br label %HUF_decompress1X2_usingDTable_internal.exit

HUF_decompress1X2_usingDTable_internal.exit:      ; preds = %cond.false, %cleanup.cont.i
  %retval.1.i = phi i64 [ %call2.i, %cond.false ], [ %.dstSize.i, %cleanup.cont.i ]
  call void @llvm.lifetime.end(i64 32, i8* %2) #11
  br label %cond.end

cond.end:                                         ; preds = %HUF_decompress1X2_usingDTable_internal.exit, %HUF_decompress1X4_usingDTable_internal.exit
  %cond = phi i64 [ %retval.1.i19, %HUF_decompress1X4_usingDTable_internal.exit ], [ %retval.1.i, %HUF_decompress1X2_usingDTable_internal.exit ]
  ret i64 %cond
}

; Function Attrs: nounwind ssp uwtable
define i64 @HUF_decompress4X_usingDTable(i8* %dst, i64 %maxDstSize, i8* %cSrc, i64 %cSrcSize, i32* %DTable) #0 {
entry:
  %DTable.val = load i32, i32* %DTable, align 1
  %conv = and i32 %DTable.val, 65280
  %tobool = icmp eq i32 %conv, 0
  br i1 %tobool, label %cond.false, label %cond.true

cond.true:                                        ; preds = %entry
  %call1 = tail call fastcc i64 @HUF_decompress4X4_usingDTable_internal(i8* %dst, i64 %maxDstSize, i8* %cSrc, i64 %cSrcSize, i32* nonnull %DTable)
  br label %cond.end

cond.false:                                       ; preds = %entry
  %call2 = tail call fastcc i64 @HUF_decompress4X2_usingDTable_internal(i8* %dst, i64 %maxDstSize, i8* %cSrc, i64 %cSrcSize, i32* nonnull %DTable)
  br label %cond.end

cond.end:                                         ; preds = %cond.false, %cond.true
  %cond = phi i64 [ %call1, %cond.true ], [ %call2, %cond.false ]
  ret i64 %cond
}

; Function Attrs: norecurse nounwind readnone ssp uwtable
define i32 @HUF_selectDecoder(i64 %dstSize, i64 %cSrcSize) #5 {
entry:
  %mul = shl i64 %cSrcSize, 4
  %div = udiv i64 %mul, %dstSize
  %shr = lshr i64 %dstSize, 8
  %conv1 = trunc i64 %shr to i32
  %idxprom = and i64 %div, 4294967295
  %tableTime = getelementptr inbounds [16 x [3 x %struct.algo_time_t]], [16 x [3 x %struct.algo_time_t]]* @algoTime, i64 0, i64 %idxprom, i64 0, i32 0
  %0 = load i32, i32* %tableTime, align 8, !tbaa !36
  %decode256Time = getelementptr inbounds [16 x [3 x %struct.algo_time_t]], [16 x [3 x %struct.algo_time_t]]* @algoTime, i64 0, i64 %idxprom, i64 0, i32 1
  %1 = load i32, i32* %decode256Time, align 4, !tbaa !38
  %mul6 = mul i32 %1, %conv1
  %add = add i32 %mul6, %0
  %tableTime10 = getelementptr inbounds [16 x [3 x %struct.algo_time_t]], [16 x [3 x %struct.algo_time_t]]* @algoTime, i64 0, i64 %idxprom, i64 1, i32 0
  %2 = load i32, i32* %tableTime10, align 8, !tbaa !36
  %decode256Time14 = getelementptr inbounds [16 x [3 x %struct.algo_time_t]], [16 x [3 x %struct.algo_time_t]]* @algoTime, i64 0, i64 %idxprom, i64 1, i32 1
  %3 = load i32, i32* %decode256Time14, align 4, !tbaa !38
  %mul15 = mul i32 %3, %conv1
  %add16 = add i32 %mul15, %2
  %shr17 = lshr i32 %add16, 3
  %add18 = add i32 %shr17, %add16
  %cmp = icmp ult i32 %add18, %add
  %conv19 = zext i1 %cmp to i32
  ret i32 %conv19
}

; Function Attrs: nounwind ssp uwtable
define i64 @HUF_decompress(i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize) #0 {
entry:
  %cmp = icmp eq i64 %dstSize, 0
  br i1 %cmp, label %return, label %if.end

if.end:                                           ; preds = %entry
  %cmp1 = icmp ugt i64 %cSrcSize, %dstSize
  br i1 %cmp1, label %return, label %if.end3

if.end3:                                          ; preds = %if.end
  %cmp4 = icmp eq i64 %cSrcSize, %dstSize
  br i1 %cmp4, label %if.then5, label %if.end6

if.then5:                                         ; preds = %if.end3
  %0 = tail call i64 @llvm.objectsize.i64.p0i8(i8* %dst, i1 false)
  %call = tail call i8* @__memcpy_chk(i8* %dst, i8* %cSrc, i64 %dstSize, i64 %0) #11
  br label %return

if.end6:                                          ; preds = %if.end3
  %cmp7 = icmp eq i64 %cSrcSize, 1
  br i1 %cmp7, label %if.then8, label %if.end10

if.then8:                                         ; preds = %if.end6
  %1 = load i8, i8* %cSrc, align 1, !tbaa !11
  %conv = zext i8 %1 to i32
  %2 = tail call i64 @llvm.objectsize.i64.p0i8(i8* %dst, i1 false)
  %call9 = tail call i8* @__memset_chk(i8* %dst, i32 %conv, i64 %dstSize, i64 %2) #11
  br label %return

if.end10:                                         ; preds = %if.end6
  %mul.i = shl i64 %cSrcSize, 4
  %div.i = udiv i64 %mul.i, %dstSize
  %shr.i = lshr i64 %dstSize, 8
  %conv1.i = trunc i64 %shr.i to i32
  %idxprom.i = and i64 %div.i, 4294967295
  %tableTime.i = getelementptr inbounds [16 x [3 x %struct.algo_time_t]], [16 x [3 x %struct.algo_time_t]]* @algoTime, i64 0, i64 %idxprom.i, i64 0, i32 0
  %3 = load i32, i32* %tableTime.i, align 8, !tbaa !36
  %decode256Time.i = getelementptr inbounds [16 x [3 x %struct.algo_time_t]], [16 x [3 x %struct.algo_time_t]]* @algoTime, i64 0, i64 %idxprom.i, i64 0, i32 1
  %4 = load i32, i32* %decode256Time.i, align 4, !tbaa !38
  %mul6.i = mul i32 %4, %conv1.i
  %add.i = add i32 %mul6.i, %3
  %tableTime10.i = getelementptr inbounds [16 x [3 x %struct.algo_time_t]], [16 x [3 x %struct.algo_time_t]]* @algoTime, i64 0, i64 %idxprom.i, i64 1, i32 0
  %5 = load i32, i32* %tableTime10.i, align 8, !tbaa !36
  %decode256Time14.i = getelementptr inbounds [16 x [3 x %struct.algo_time_t]], [16 x [3 x %struct.algo_time_t]]* @algoTime, i64 0, i64 %idxprom.i, i64 1, i32 1
  %6 = load i32, i32* %decode256Time14.i, align 4, !tbaa !38
  %mul15.i = mul i32 %6, %conv1.i
  %add16.i = add i32 %mul15.i, %5
  %shr17.i = lshr i32 %add16.i, 3
  %add18.i = add i32 %shr17.i, %add16.i
  %cmp.i = icmp ult i32 %add18.i, %add.i
  %idxprom = zext i1 %cmp.i to i64
  %arrayidx = getelementptr inbounds [2 x i64 (i8*, i64, i8*, i64)*], [2 x i64 (i8*, i64, i8*, i64)*]* @HUF_decompress.decompress, i64 0, i64 %idxprom
  %7 = load i64 (i8*, i64, i8*, i64)*, i64 (i8*, i64, i8*, i64)** %arrayidx, align 8, !tbaa !39
  %call12 = tail call i64 %7(i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize) #11
  br label %return

return:                                           ; preds = %if.end, %entry, %if.end10, %if.then8, %if.then5
  %retval.0 = phi i64 [ %dstSize, %if.then5 ], [ %dstSize, %if.then8 ], [ %call12, %if.end10 ], [ -12, %entry ], [ -14, %if.end ]
  ret i64 %retval.0
}

; Function Attrs: nounwind
declare i8* @__memset_chk(i8*, i32, i64, i64) #3

; Function Attrs: nounwind ssp uwtable
define i64 @HUF_decompress4X_DCtx(i32* %dctx, i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize) #0 {
entry:
  %cmp = icmp eq i64 %dstSize, 0
  br i1 %cmp, label %return, label %if.end

if.end:                                           ; preds = %entry
  %cmp1 = icmp ugt i64 %cSrcSize, %dstSize
  br i1 %cmp1, label %return, label %if.end3

if.end3:                                          ; preds = %if.end
  %cmp4 = icmp eq i64 %cSrcSize, %dstSize
  br i1 %cmp4, label %if.then5, label %if.end6

if.then5:                                         ; preds = %if.end3
  %0 = tail call i64 @llvm.objectsize.i64.p0i8(i8* %dst, i1 false)
  %call = tail call i8* @__memcpy_chk(i8* %dst, i8* %cSrc, i64 %dstSize, i64 %0) #11
  br label %return

if.end6:                                          ; preds = %if.end3
  %cmp7 = icmp eq i64 %cSrcSize, 1
  br i1 %cmp7, label %if.then8, label %if.end10

if.then8:                                         ; preds = %if.end6
  %1 = load i8, i8* %cSrc, align 1, !tbaa !11
  %conv = zext i8 %1 to i32
  %2 = tail call i64 @llvm.objectsize.i64.p0i8(i8* %dst, i1 false)
  %call9 = tail call i8* @__memset_chk(i8* %dst, i32 %conv, i64 %dstSize, i64 %2) #11
  br label %return

if.end10:                                         ; preds = %if.end6
  %mul.i = shl i64 %cSrcSize, 4
  %div.i = udiv i64 %mul.i, %dstSize
  %shr.i = lshr i64 %dstSize, 8
  %conv1.i = trunc i64 %shr.i to i32
  %idxprom.i = and i64 %div.i, 4294967295
  %tableTime.i = getelementptr inbounds [16 x [3 x %struct.algo_time_t]], [16 x [3 x %struct.algo_time_t]]* @algoTime, i64 0, i64 %idxprom.i, i64 0, i32 0
  %3 = load i32, i32* %tableTime.i, align 8, !tbaa !36
  %decode256Time.i = getelementptr inbounds [16 x [3 x %struct.algo_time_t]], [16 x [3 x %struct.algo_time_t]]* @algoTime, i64 0, i64 %idxprom.i, i64 0, i32 1
  %4 = load i32, i32* %decode256Time.i, align 4, !tbaa !38
  %mul6.i = mul i32 %4, %conv1.i
  %add.i = add i32 %mul6.i, %3
  %tableTime10.i = getelementptr inbounds [16 x [3 x %struct.algo_time_t]], [16 x [3 x %struct.algo_time_t]]* @algoTime, i64 0, i64 %idxprom.i, i64 1, i32 0
  %5 = load i32, i32* %tableTime10.i, align 8, !tbaa !36
  %decode256Time14.i = getelementptr inbounds [16 x [3 x %struct.algo_time_t]], [16 x [3 x %struct.algo_time_t]]* @algoTime, i64 0, i64 %idxprom.i, i64 1, i32 1
  %6 = load i32, i32* %decode256Time14.i, align 4, !tbaa !38
  %mul15.i = mul i32 %6, %conv1.i
  %add16.i = add i32 %mul15.i, %5
  %shr17.i = lshr i32 %add16.i, 3
  %add18.i = add i32 %shr17.i, %add16.i
  %cmp.i = icmp ult i32 %add18.i, %add.i
  br i1 %cmp.i, label %cond.true, label %cond.false

cond.true:                                        ; preds = %if.end10
  %call.i34 = tail call i64 @HUF_readDTableX4(i32* %dctx, i8* %cSrc, i64 %cSrcSize) #11
  %call1.i35 = tail call i32 @HUF_isError(i64 %call.i34) #11
  %tobool.i36 = icmp eq i32 %call1.i35, 0
  br i1 %tobool.i36, label %if.end.i38, label %return

if.end.i38:                                       ; preds = %cond.true
  %cmp.i37 = icmp ult i64 %call.i34, %cSrcSize
  br i1 %cmp.i37, label %if.end3.i42, label %return

if.end3.i42:                                      ; preds = %if.end.i38
  %add.ptr.i39 = getelementptr inbounds i8, i8* %cSrc, i64 %call.i34
  %sub.i40 = sub i64 %cSrcSize, %call.i34
  %call4.i41 = tail call fastcc i64 @HUF_decompress4X4_usingDTable_internal(i8* %dst, i64 %dstSize, i8* %add.ptr.i39, i64 %sub.i40, i32* %dctx) #11
  br label %return

cond.false:                                       ; preds = %if.end10
  %call.i = tail call i64 @HUF_readDTableX2(i32* %dctx, i8* %cSrc, i64 %cSrcSize) #11
  %call1.i = tail call i32 @HUF_isError(i64 %call.i) #11
  %tobool.i = icmp eq i32 %call1.i, 0
  br i1 %tobool.i, label %if.end.i, label %return

if.end.i:                                         ; preds = %cond.false
  %cmp.i33 = icmp ult i64 %call.i, %cSrcSize
  br i1 %cmp.i33, label %if.end3.i, label %return

if.end3.i:                                        ; preds = %if.end.i
  %add.ptr.i = getelementptr inbounds i8, i8* %cSrc, i64 %call.i
  %sub.i = sub i64 %cSrcSize, %call.i
  %call4.i = tail call fastcc i64 @HUF_decompress4X2_usingDTable_internal(i8* %dst, i64 %dstSize, i8* %add.ptr.i, i64 %sub.i, i32* %dctx) #11
  br label %return

return:                                           ; preds = %if.end3.i, %if.end.i, %cond.false, %if.end3.i42, %if.end.i38, %cond.true, %if.end, %entry, %if.then8, %if.then5
  %retval.0 = phi i64 [ %dstSize, %if.then5 ], [ %dstSize, %if.then8 ], [ -12, %entry ], [ -14, %if.end ], [ %call4.i41, %if.end3.i42 ], [ %call.i34, %cond.true ], [ -13, %if.end.i38 ], [ %call4.i, %if.end3.i ], [ %call.i, %cond.false ], [ -13, %if.end.i ]
  ret i64 %retval.0
}

; Function Attrs: nounwind ssp uwtable
define i64 @HUF_decompress4X_hufOnly(i32* %dctx, i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize) #0 {
entry:
  %cmp = icmp eq i64 %dstSize, 0
  br i1 %cmp, label %return, label %if.end

if.end:                                           ; preds = %entry
  %cmp1 = icmp uge i64 %cSrcSize, %dstSize
  %cmp2 = icmp ult i64 %cSrcSize, 2
  %or.cond = or i1 %cmp1, %cmp2
  br i1 %or.cond, label %return, label %if.end4

if.end4:                                          ; preds = %if.end
  %mul.i = shl i64 %cSrcSize, 4
  %div.i = udiv i64 %mul.i, %dstSize
  %shr.i = lshr i64 %dstSize, 8
  %conv1.i = trunc i64 %shr.i to i32
  %idxprom.i = and i64 %div.i, 4294967295
  %tableTime.i = getelementptr inbounds [16 x [3 x %struct.algo_time_t]], [16 x [3 x %struct.algo_time_t]]* @algoTime, i64 0, i64 %idxprom.i, i64 0, i32 0
  %0 = load i32, i32* %tableTime.i, align 8, !tbaa !36
  %decode256Time.i = getelementptr inbounds [16 x [3 x %struct.algo_time_t]], [16 x [3 x %struct.algo_time_t]]* @algoTime, i64 0, i64 %idxprom.i, i64 0, i32 1
  %1 = load i32, i32* %decode256Time.i, align 4, !tbaa !38
  %mul6.i = mul i32 %1, %conv1.i
  %add.i = add i32 %mul6.i, %0
  %tableTime10.i = getelementptr inbounds [16 x [3 x %struct.algo_time_t]], [16 x [3 x %struct.algo_time_t]]* @algoTime, i64 0, i64 %idxprom.i, i64 1, i32 0
  %2 = load i32, i32* %tableTime10.i, align 8, !tbaa !36
  %decode256Time14.i = getelementptr inbounds [16 x [3 x %struct.algo_time_t]], [16 x [3 x %struct.algo_time_t]]* @algoTime, i64 0, i64 %idxprom.i, i64 1, i32 1
  %3 = load i32, i32* %decode256Time14.i, align 4, !tbaa !38
  %mul15.i = mul i32 %3, %conv1.i
  %add16.i = add i32 %mul15.i, %2
  %shr17.i = lshr i32 %add16.i, 3
  %add18.i = add i32 %shr17.i, %add16.i
  %cmp.i = icmp ult i32 %add18.i, %add.i
  br i1 %cmp.i, label %cond.true, label %cond.false

cond.true:                                        ; preds = %if.end4
  %call.i15 = tail call i64 @HUF_readDTableX4(i32* %dctx, i8* %cSrc, i64 %cSrcSize) #11
  %call1.i16 = tail call i32 @HUF_isError(i64 %call.i15) #11
  %tobool.i17 = icmp eq i32 %call1.i16, 0
  br i1 %tobool.i17, label %if.end.i19, label %return

if.end.i19:                                       ; preds = %cond.true
  %cmp.i18 = icmp ult i64 %call.i15, %cSrcSize
  br i1 %cmp.i18, label %if.end3.i23, label %return

if.end3.i23:                                      ; preds = %if.end.i19
  %add.ptr.i20 = getelementptr inbounds i8, i8* %cSrc, i64 %call.i15
  %sub.i21 = sub i64 %cSrcSize, %call.i15
  %call4.i22 = tail call fastcc i64 @HUF_decompress4X4_usingDTable_internal(i8* %dst, i64 %dstSize, i8* %add.ptr.i20, i64 %sub.i21, i32* %dctx) #11
  br label %return

cond.false:                                       ; preds = %if.end4
  %call.i = tail call i64 @HUF_readDTableX2(i32* %dctx, i8* %cSrc, i64 %cSrcSize) #11
  %call1.i = tail call i32 @HUF_isError(i64 %call.i) #11
  %tobool.i = icmp eq i32 %call1.i, 0
  br i1 %tobool.i, label %if.end.i, label %return

if.end.i:                                         ; preds = %cond.false
  %cmp.i14 = icmp ult i64 %call.i, %cSrcSize
  br i1 %cmp.i14, label %if.end3.i, label %return

if.end3.i:                                        ; preds = %if.end.i
  %add.ptr.i = getelementptr inbounds i8, i8* %cSrc, i64 %call.i
  %sub.i = sub i64 %cSrcSize, %call.i
  %call4.i = tail call fastcc i64 @HUF_decompress4X2_usingDTable_internal(i8* %dst, i64 %dstSize, i8* %add.ptr.i, i64 %sub.i, i32* %dctx) #11
  br label %return

return:                                           ; preds = %if.end3.i, %if.end.i, %cond.false, %if.end3.i23, %if.end.i19, %cond.true, %if.end, %entry
  %retval.0 = phi i64 [ -12, %entry ], [ -14, %if.end ], [ %call4.i22, %if.end3.i23 ], [ %call.i15, %cond.true ], [ -13, %if.end.i19 ], [ %call4.i, %if.end3.i ], [ %call.i, %cond.false ], [ -13, %if.end.i ]
  ret i64 %retval.0
}

; Function Attrs: nounwind ssp uwtable
define i64 @HUF_decompress1X_DCtx(i32* %dctx, i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize) #0 {
entry:
  %cmp = icmp eq i64 %dstSize, 0
  br i1 %cmp, label %return, label %if.end

if.end:                                           ; preds = %entry
  %cmp1 = icmp ugt i64 %cSrcSize, %dstSize
  br i1 %cmp1, label %return, label %if.end3

if.end3:                                          ; preds = %if.end
  %cmp4 = icmp eq i64 %cSrcSize, %dstSize
  br i1 %cmp4, label %if.then5, label %if.end6

if.then5:                                         ; preds = %if.end3
  %0 = tail call i64 @llvm.objectsize.i64.p0i8(i8* %dst, i1 false)
  %call = tail call i8* @__memcpy_chk(i8* %dst, i8* %cSrc, i64 %dstSize, i64 %0) #11
  br label %return

if.end6:                                          ; preds = %if.end3
  %cmp7 = icmp eq i64 %cSrcSize, 1
  br i1 %cmp7, label %if.then8, label %if.end10

if.then8:                                         ; preds = %if.end6
  %1 = load i8, i8* %cSrc, align 1, !tbaa !11
  %conv = zext i8 %1 to i32
  %2 = tail call i64 @llvm.objectsize.i64.p0i8(i8* %dst, i1 false)
  %call9 = tail call i8* @__memset_chk(i8* %dst, i32 %conv, i64 %dstSize, i64 %2) #11
  br label %return

if.end10:                                         ; preds = %if.end6
  %mul.i = shl i64 %cSrcSize, 4
  %div.i = udiv i64 %mul.i, %dstSize
  %shr.i = lshr i64 %dstSize, 8
  %conv1.i = trunc i64 %shr.i to i32
  %idxprom.i = and i64 %div.i, 4294967295
  %tableTime.i = getelementptr inbounds [16 x [3 x %struct.algo_time_t]], [16 x [3 x %struct.algo_time_t]]* @algoTime, i64 0, i64 %idxprom.i, i64 0, i32 0
  %3 = load i32, i32* %tableTime.i, align 8, !tbaa !36
  %decode256Time.i = getelementptr inbounds [16 x [3 x %struct.algo_time_t]], [16 x [3 x %struct.algo_time_t]]* @algoTime, i64 0, i64 %idxprom.i, i64 0, i32 1
  %4 = load i32, i32* %decode256Time.i, align 4, !tbaa !38
  %mul6.i = mul i32 %4, %conv1.i
  %add.i = add i32 %mul6.i, %3
  %tableTime10.i = getelementptr inbounds [16 x [3 x %struct.algo_time_t]], [16 x [3 x %struct.algo_time_t]]* @algoTime, i64 0, i64 %idxprom.i, i64 1, i32 0
  %5 = load i32, i32* %tableTime10.i, align 8, !tbaa !36
  %decode256Time14.i = getelementptr inbounds [16 x [3 x %struct.algo_time_t]], [16 x [3 x %struct.algo_time_t]]* @algoTime, i64 0, i64 %idxprom.i, i64 1, i32 1
  %6 = load i32, i32* %decode256Time14.i, align 4, !tbaa !38
  %mul15.i = mul i32 %6, %conv1.i
  %add16.i = add i32 %mul15.i, %5
  %shr17.i = lshr i32 %add16.i, 3
  %add18.i = add i32 %shr17.i, %add16.i
  %cmp.i = icmp ult i32 %add18.i, %add.i
  br i1 %cmp.i, label %cond.true, label %cond.false

cond.true:                                        ; preds = %if.end10
  %call12 = tail call i64 @HUF_decompress1X4_DCtx(i32* %dctx, i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize)
  br label %return

cond.false:                                       ; preds = %if.end10
  %call13 = tail call i64 @HUF_decompress1X2_DCtx(i32* %dctx, i8* %dst, i64 %dstSize, i8* %cSrc, i64 %cSrcSize)
  br label %return

return:                                           ; preds = %cond.true, %cond.false, %if.end, %entry, %if.then8, %if.then5
  %retval.0 = phi i64 [ %dstSize, %if.then5 ], [ %dstSize, %if.then8 ], [ -12, %entry ], [ -14, %if.end ], [ %call12, %cond.true ], [ %call13, %cond.false ]
  ret i64 %retval.0
}

; Function Attrs: inlinehint nounwind ssp uwtable
define internal fastcc i64 @BIT_initDStream(%struct.BIT_DStream_t* %bitD, i8* %srcBuffer, i64 %srcSize) unnamed_addr #6 {
entry:
  %cmp = icmp eq i64 %srcSize, 0
  br i1 %cmp, label %if.then, label %if.end

if.then:                                          ; preds = %entry
  %0 = bitcast %struct.BIT_DStream_t* %bitD to i8*
  %1 = tail call i64 @llvm.objectsize.i64.p0i8(i8* %0, i1 false)
  %call = tail call i8* @__memset_chk(i8* %0, i32 0, i64 32, i64 %1) #11
  br label %return

if.end:                                           ; preds = %entry
  %cmp1 = icmp ugt i64 %srcSize, 7
  %start = getelementptr inbounds %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD, i64 0, i32 3
  store i8* %srcBuffer, i8** %start, align 8, !tbaa !18
  br i1 %cmp1, label %if.then2, label %if.else

if.then2:                                         ; preds = %if.end
  %add.ptr = getelementptr inbounds i8, i8* %srcBuffer, i64 %srcSize
  %add.ptr3 = getelementptr inbounds i8, i8* %add.ptr, i64 -8
  %ptr = getelementptr inbounds %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD, i64 0, i32 2
  store i8* %add.ptr3, i8** %ptr, align 8, !tbaa !17
  %val.0..sroa_cast1.i.i.i = bitcast i8* %add.ptr3 to i64*
  %val.0.copyload.i.i.i = load i64, i64* %val.0..sroa_cast1.i.i.i, align 1
  %bitContainer = getelementptr inbounds %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD, i64 0, i32 0
  store i64 %val.0.copyload.i.i.i, i64* %bitContainer, align 8, !tbaa !40
  %sub = add i64 %srcSize, -1
  %arrayidx = getelementptr inbounds i8, i8* %srcBuffer, i64 %sub
  %2 = load i8, i8* %arrayidx, align 1, !tbaa !11
  %tobool = icmp eq i8 %2, 0
  br i1 %tobool, label %cond.end.thread123, label %cond.end.thread

cond.end.thread123:                               ; preds = %if.then2
  %bitsConsumed124 = getelementptr inbounds %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD, i64 0, i32 1
  store i32 0, i32* %bitsConsumed124, align 8, !tbaa !13
  br label %return

cond.end.thread:                                  ; preds = %if.then2
  %conv = zext i8 %2 to i32
  %3 = tail call i32 @llvm.ctlz.i32(i32 %conv, i1 true) #11
  %sub.i120 = xor i32 %3, 31
  %sub8 = sub nsw i32 8, %sub.i120
  %bitsConsumed121 = getelementptr inbounds %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD, i64 0, i32 1
  store i32 %sub8, i32* %bitsConsumed121, align 8, !tbaa !13
  br label %return

if.else:                                          ; preds = %if.end
  %ptr16 = getelementptr inbounds %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD, i64 0, i32 2
  store i8* %srcBuffer, i8** %ptr16, align 8, !tbaa !17
  %4 = load i8, i8* %srcBuffer, align 1, !tbaa !11
  %conv18 = zext i8 %4 to i64
  %bitContainer19 = getelementptr inbounds %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD, i64 0, i32 0
  store i64 %conv18, i64* %bitContainer19, align 8, !tbaa !40
  switch i64 %srcSize, label %sw.epilog [
    i64 7, label %sw.bb
    i64 6, label %sw.bb23
    i64 5, label %sw.bb29
    i64 4, label %sw.bb35
    i64 3, label %sw.bb41
    i64 2, label %sw.bb47
  ]

sw.bb:                                            ; preds = %if.else
  %arrayidx20 = getelementptr inbounds i8, i8* %srcBuffer, i64 6
  %5 = load i8, i8* %arrayidx20, align 1, !tbaa !11
  %conv21 = zext i8 %5 to i64
  %shl = shl nuw nsw i64 %conv21, 48
  %add = or i64 %shl, %conv18
  store i64 %add, i64* %bitContainer19, align 8, !tbaa !40
  br label %sw.bb23

sw.bb23:                                          ; preds = %if.else, %sw.bb
  %6 = phi i64 [ %conv18, %if.else ], [ %add, %sw.bb ]
  %arrayidx24 = getelementptr inbounds i8, i8* %srcBuffer, i64 5
  %7 = load i8, i8* %arrayidx24, align 1, !tbaa !11
  %conv25 = zext i8 %7 to i64
  %shl26 = shl nuw nsw i64 %conv25, 40
  %add28 = add nsw i64 %shl26, %6
  store i64 %add28, i64* %bitContainer19, align 8, !tbaa !40
  br label %sw.bb29

sw.bb29:                                          ; preds = %if.else, %sw.bb23
  %8 = phi i64 [ %conv18, %if.else ], [ %add28, %sw.bb23 ]
  %arrayidx30 = getelementptr inbounds i8, i8* %srcBuffer, i64 4
  %9 = load i8, i8* %arrayidx30, align 1, !tbaa !11
  %conv31 = zext i8 %9 to i64
  %shl32 = shl nuw nsw i64 %conv31, 32
  %add34 = add nsw i64 %shl32, %8
  store i64 %add34, i64* %bitContainer19, align 8, !tbaa !40
  br label %sw.bb35

sw.bb35:                                          ; preds = %if.else, %sw.bb29
  %10 = phi i64 [ %conv18, %if.else ], [ %add34, %sw.bb29 ]
  %arrayidx36 = getelementptr inbounds i8, i8* %srcBuffer, i64 3
  %11 = load i8, i8* %arrayidx36, align 1, !tbaa !11
  %conv37 = zext i8 %11 to i64
  %shl38 = shl nuw nsw i64 %conv37, 24
  %add40 = add i64 %shl38, %10
  store i64 %add40, i64* %bitContainer19, align 8, !tbaa !40
  br label %sw.bb41

sw.bb41:                                          ; preds = %if.else, %sw.bb35
  %12 = phi i64 [ %conv18, %if.else ], [ %add40, %sw.bb35 ]
  %arrayidx42 = getelementptr inbounds i8, i8* %srcBuffer, i64 2
  %13 = load i8, i8* %arrayidx42, align 1, !tbaa !11
  %conv43 = zext i8 %13 to i64
  %shl44 = shl nuw nsw i64 %conv43, 16
  %add46 = add i64 %shl44, %12
  store i64 %add46, i64* %bitContainer19, align 8, !tbaa !40
  br label %sw.bb47

sw.bb47:                                          ; preds = %if.else, %sw.bb41
  %14 = phi i64 [ %conv18, %if.else ], [ %add46, %sw.bb41 ]
  %arrayidx48 = getelementptr inbounds i8, i8* %srcBuffer, i64 1
  %15 = load i8, i8* %arrayidx48, align 1, !tbaa !11
  %conv49 = zext i8 %15 to i64
  %shl50 = shl nuw nsw i64 %conv49, 8
  %add52 = add i64 %shl50, %14
  store i64 %add52, i64* %bitContainer19, align 8, !tbaa !40
  br label %sw.epilog

sw.epilog:                                        ; preds = %sw.bb47, %if.else
  %sub54 = add i64 %srcSize, -1
  %arrayidx55 = getelementptr inbounds i8, i8* %srcBuffer, i64 %sub54
  %16 = load i8, i8* %arrayidx55, align 1, !tbaa !11
  %tobool57 = icmp eq i8 %16, 0
  br i1 %tobool57, label %cond.end63.thread130, label %cleanup.cont73

cond.end63.thread130:                             ; preds = %sw.epilog
  %bitsConsumed65131 = getelementptr inbounds %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD, i64 0, i32 1
  store i32 0, i32* %bitsConsumed65131, align 8, !tbaa !13
  br label %return

cleanup.cont73:                                   ; preds = %sw.epilog
  %conv56 = zext i8 %16 to i32
  %17 = tail call i32 @llvm.ctlz.i32(i32 %conv56, i1 true) #11
  %sub.i = xor i32 %17, 31
  %bitsConsumed65126 = getelementptr inbounds %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD, i64 0, i32 1
  %sub74 = sub i64 8, %srcSize
  %conv75 = trunc i64 %sub74 to i32
  %mul = shl i32 %conv75, 3
  %sub61 = add i32 %mul, 8
  %add77 = sub i32 %sub61, %sub.i
  store i32 %add77, i32* %bitsConsumed65126, align 8, !tbaa !13
  br label %return

return:                                           ; preds = %cleanup.cont73, %cond.end.thread, %cond.end63.thread130, %cond.end.thread123, %if.then
  %retval.2 = phi i64 [ -13, %if.then ], [ -1, %cond.end.thread123 ], [ -1, %cond.end63.thread130 ], [ %srcSize, %cond.end.thread ], [ %srcSize, %cleanup.cont73 ]
  ret i64 %retval.2
}

; Function Attrs: inlinehint norecurse nounwind ssp uwtable
define internal fastcc void @HUF_decodeStreamX2(i8* %p, %struct.BIT_DStream_t* nocapture %bitDPtr, i8* readnone %pEnd, %struct.HUF_DEltX2* nocapture readonly %dt, i32 %dtLog) unnamed_addr #7 {
entry:
  %call9 = tail call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* %bitDPtr)
  %cmp10 = icmp ne i32 %call9, 0
  %add.ptr = getelementptr inbounds i8, i8* %pEnd, i64 -4
  %cmp111 = icmp ult i8* %add.ptr, %p
  %or.cond12 = or i1 %cmp111, %cmp10
  br i1 %or.cond12, label %while.cond18.preheader, label %while.body.preheader

while.body.preheader:                             ; preds = %entry
  br label %while.body

while.cond18.preheader.loopexit:                  ; preds = %while.body
  %incdec.ptr17.lcssa = phi i8* [ %incdec.ptr17, %while.body ]
  br label %while.cond18.preheader

while.cond18.preheader:                           ; preds = %while.cond18.preheader.loopexit, %entry
  %p.addr.0.lcssa = phi i8* [ %p, %entry ], [ %incdec.ptr17.lcssa, %while.cond18.preheader.loopexit ]
  %call194 = tail call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* %bitDPtr)
  %cmp205 = icmp eq i32 %call194, 0
  %cmp226 = icmp ult i8* %p.addr.0.lcssa, %pEnd
  %or.cond17 = and i1 %cmp226, %cmp205
  br i1 %or.cond17, label %while.body24.preheader, label %while.cond28.preheader

while.body24.preheader:                           ; preds = %while.cond18.preheader
  br label %while.body24

while.body:                                       ; preds = %while.body.preheader, %while.body
  %p.addr.013 = phi i8* [ %incdec.ptr17, %while.body ], [ %p, %while.body.preheader ]
  %call3 = tail call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* %bitDPtr, %struct.HUF_DEltX2* %dt, i32 %dtLog)
  %incdec.ptr = getelementptr inbounds i8, i8* %p.addr.013, i64 1
  store i8 %call3, i8* %p.addr.013, align 1, !tbaa !11
  %call7 = tail call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* %bitDPtr, %struct.HUF_DEltX2* %dt, i32 %dtLog)
  %incdec.ptr8 = getelementptr inbounds i8, i8* %p.addr.013, i64 2
  store i8 %call7, i8* %incdec.ptr, align 1, !tbaa !11
  %call13 = tail call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* %bitDPtr, %struct.HUF_DEltX2* %dt, i32 %dtLog)
  %incdec.ptr14 = getelementptr inbounds i8, i8* %p.addr.013, i64 3
  store i8 %call13, i8* %incdec.ptr8, align 1, !tbaa !11
  %call16 = tail call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* %bitDPtr, %struct.HUF_DEltX2* %dt, i32 %dtLog)
  %incdec.ptr17 = getelementptr inbounds i8, i8* %p.addr.013, i64 4
  store i8 %call16, i8* %incdec.ptr14, align 1, !tbaa !11
  %call = tail call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* %bitDPtr)
  %cmp = icmp ne i32 %call, 0
  %cmp1 = icmp ugt i8* %incdec.ptr17, %add.ptr
  %or.cond = or i1 %cmp1, %cmp
  br i1 %or.cond, label %while.cond18.preheader.loopexit, label %while.body

while.cond28.preheader.loopexit:                  ; preds = %while.body24
  %incdec.ptr26.lcssa = phi i8* [ %incdec.ptr26, %while.body24 ]
  br label %while.cond28.preheader

while.cond28.preheader:                           ; preds = %while.cond28.preheader.loopexit, %while.cond18.preheader
  %p.addr.3.lcssa = phi i8* [ %p.addr.0.lcssa, %while.cond18.preheader ], [ %incdec.ptr26.lcssa, %while.cond28.preheader.loopexit ]
  %cmp292 = icmp ult i8* %p.addr.3.lcssa, %pEnd
  br i1 %cmp292, label %while.body30.preheader, label %while.end33

while.body30.preheader:                           ; preds = %while.cond28.preheader
  br label %while.body30

while.body24:                                     ; preds = %while.body24.preheader, %while.body24
  %p.addr.38 = phi i8* [ %incdec.ptr26, %while.body24 ], [ %p.addr.0.lcssa, %while.body24.preheader ]
  %call25 = tail call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* %bitDPtr, %struct.HUF_DEltX2* %dt, i32 %dtLog)
  %incdec.ptr26 = getelementptr inbounds i8, i8* %p.addr.38, i64 1
  store i8 %call25, i8* %p.addr.38, align 1, !tbaa !11
  %call19 = tail call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* %bitDPtr)
  %cmp20 = icmp eq i32 %call19, 0
  %cmp22 = icmp ult i8* %incdec.ptr26, %pEnd
  %or.cond1 = and i1 %cmp22, %cmp20
  br i1 %or.cond1, label %while.body24, label %while.cond28.preheader.loopexit

while.body30:                                     ; preds = %while.body30.preheader, %while.body30
  %p.addr.43 = phi i8* [ %incdec.ptr32, %while.body30 ], [ %p.addr.3.lcssa, %while.body30.preheader ]
  %call31 = tail call fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* %bitDPtr, %struct.HUF_DEltX2* %dt, i32 %dtLog)
  %incdec.ptr32 = getelementptr inbounds i8, i8* %p.addr.43, i64 1
  store i8 %call31, i8* %p.addr.43, align 1, !tbaa !11
  %exitcond = icmp eq i8* %incdec.ptr32, %pEnd
  br i1 %exitcond, label %while.end33.loopexit, label %while.body30

while.end33.loopexit:                             ; preds = %while.body30
  br label %while.end33

while.end33:                                      ; preds = %while.end33.loopexit, %while.cond28.preheader
  ret void
}

; Function Attrs: nounwind readnone
declare i32 @llvm.ctlz.i32(i32, i1) #4

; Function Attrs: inlinehint noinline norecurse nounwind ssp uwtable
define internal fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* nocapture %bitD) unnamed_addr #8 {
entry:
  %bitsConsumed = getelementptr inbounds %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD, i64 0, i32 1
  %0 = load i32, i32* %bitsConsumed, align 8, !tbaa !13
  %cmp = icmp ugt i32 %0, 64
  br i1 %cmp, label %return, label %if.end

if.end:                                           ; preds = %entry
  %ptr = getelementptr inbounds %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD, i64 0, i32 2
  %1 = load i8*, i8** %ptr, align 8, !tbaa !17
  %start = getelementptr inbounds %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD, i64 0, i32 3
  %2 = load i8*, i8** %start, align 8, !tbaa !18
  %add.ptr = getelementptr inbounds i8, i8* %2, i64 8
  %cmp2 = icmp ult i8* %1, %add.ptr
  br i1 %cmp2, label %if.end10, label %if.then4

if.then4:                                         ; preds = %if.end
  %shr = lshr i32 %0, 3
  %idx.ext = zext i32 %shr to i64
  %idx.neg = sub nsw i64 0, %idx.ext
  %add.ptr7 = getelementptr inbounds i8, i8* %1, i64 %idx.neg
  store i8* %add.ptr7, i8** %ptr, align 8, !tbaa !17
  %and = and i32 %0, 7
  store i32 %and, i32* %bitsConsumed, align 8, !tbaa !13
  %val.0..sroa_cast1.i.i.i = bitcast i8* %add.ptr7 to i64*
  %val.0.copyload.i.i.i = load i64, i64* %val.0..sroa_cast1.i.i.i, align 1
  %bitContainer = getelementptr inbounds %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD, i64 0, i32 0
  store i64 %val.0.copyload.i.i.i, i64* %bitContainer, align 8, !tbaa !40
  br label %return

if.end10:                                         ; preds = %if.end
  %cmp13 = icmp eq i8* %1, %2
  br i1 %cmp13, label %if.then15, label %if.end22

if.then15:                                        ; preds = %if.end10
  %cmp18 = icmp ult i32 %0, 64
  %. = select i1 %cmp18, i32 1, i32 2
  br label %return

if.end22:                                         ; preds = %if.end10
  %shr24 = lshr i32 %0, 3
  %idx.ext26 = zext i32 %shr24 to i64
  %idx.neg27 = sub nsw i64 0, %idx.ext26
  %add.ptr28 = getelementptr inbounds i8, i8* %1, i64 %idx.neg27
  %cmp30 = icmp ult i8* %add.ptr28, %2
  %sub.ptr.lhs.cast = ptrtoint i8* %1 to i64
  %sub.ptr.rhs.cast = ptrtoint i8* %2 to i64
  %sub.ptr.sub = sub i64 %sub.ptr.lhs.cast, %sub.ptr.rhs.cast
  %conv35 = trunc i64 %sub.ptr.sub to i32
  %nbBytes.0 = select i1 %cmp30, i32 %conv35, i32 %shr24
  %result.0 = zext i1 %cmp30 to i32
  %idx.ext38 = zext i32 %nbBytes.0 to i64
  %idx.neg39 = sub nsw i64 0, %idx.ext38
  %add.ptr40 = getelementptr inbounds i8, i8* %1, i64 %idx.neg39
  store i8* %add.ptr40, i8** %ptr, align 8, !tbaa !17
  %mul = shl i32 %nbBytes.0, 3
  %sub = sub i32 %0, %mul
  store i32 %sub, i32* %bitsConsumed, align 8, !tbaa !13
  %val.0..sroa_cast1.i.i.i67 = bitcast i8* %add.ptr40 to i64*
  %val.0.copyload.i.i.i68 = load i64, i64* %val.0..sroa_cast1.i.i.i67, align 1
  %bitContainer44 = getelementptr inbounds %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitD, i64 0, i32 0
  store i64 %val.0.copyload.i.i.i68, i64* %bitContainer44, align 8, !tbaa !40
  br label %return

return:                                           ; preds = %if.then15, %entry, %if.end22, %if.then4
  %retval.0 = phi i32 [ 0, %if.then4 ], [ %result.0, %if.end22 ], [ 3, %entry ], [ %., %if.then15 ]
  ret i32 %retval.0
}

; Function Attrs: noinline norecurse nounwind ssp uwtable
define internal fastcc zeroext i8 @HUF_decodeSymbolX2(%struct.BIT_DStream_t* nocapture %Dstream, %struct.HUF_DEltX2* nocapture readonly %dt, i32 %dtLog) unnamed_addr #9 {
entry:
  %Dstream.idx = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %Dstream, i64 0, i32 0
  %Dstream.idx.val = load i64, i64* %Dstream.idx, align 8, !tbaa !40
  %Dstream.idx6 = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %Dstream, i64 0, i32 1
  %Dstream.idx6.val = load i32, i32* %Dstream.idx6, align 8, !tbaa !13
  %and.i = and i32 %Dstream.idx6.val, 63
  %sh_prom.i = zext i32 %and.i to i64
  %shl.i = shl i64 %Dstream.idx.val, %sh_prom.i
  %0 = sub i32 0, %dtLog
  %and1.i = and i32 %0, 63
  %sh_prom2.i = zext i32 %and1.i to i64
  %shr.i = lshr i64 %shl.i, %sh_prom2.i
  %byte = getelementptr inbounds %struct.HUF_DEltX2, %struct.HUF_DEltX2* %dt, i64 %shr.i, i32 0
  %1 = load i8, i8* %byte, align 1, !tbaa !20
  %nbBits = getelementptr inbounds %struct.HUF_DEltX2, %struct.HUF_DEltX2* %dt, i64 %shr.i, i32 1
  %2 = load i8, i8* %nbBits, align 1, !tbaa !22
  %conv = zext i8 %2 to i32
  %add.i = add i32 %Dstream.idx6.val, %conv
  store i32 %add.i, i32* %Dstream.idx6, align 8, !tbaa !13
  ret i8 %1
}

; Function Attrs: inlinehint nounwind ssp uwtable
define internal fastcc void @HUF_decodeStreamX4(i8* %p, %struct.BIT_DStream_t* nocapture %bitDPtr, i8* readnone %pEnd, %struct.HUF_DEltX4* %dt, i32 %dtLog) unnamed_addr #6 {
entry:
  %call12 = tail call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* %bitDPtr)
  %cmp13 = icmp eq i32 %call12, 0
  %add.ptr = getelementptr inbounds i8, i8* %pEnd, i64 -7
  %cmp114 = icmp ugt i8* %add.ptr, %p
  %and115 = and i1 %cmp114, %cmp13
  br i1 %and115, label %while.body.preheader, label %while.cond24.preheader

while.body.preheader:                             ; preds = %entry
  br label %while.body

while.cond24.preheader.loopexit:                  ; preds = %while.body
  %add.ptr23.lcssa = phi i8* [ %add.ptr23, %while.body ]
  br label %while.cond24.preheader

while.cond24.preheader:                           ; preds = %while.cond24.preheader.loopexit, %entry
  %p.addr.0.lcssa = phi i8* [ %p, %entry ], [ %add.ptr23.lcssa, %while.cond24.preheader.loopexit ]
  %call255 = tail call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* %bitDPtr)
  %cmp266 = icmp eq i32 %call255, 0
  %add.ptr28 = getelementptr inbounds i8, i8* %pEnd, i64 -2
  %cmp297 = icmp ule i8* %p.addr.0.lcssa, %add.ptr28
  %and3128 = and i1 %cmp297, %cmp266
  br i1 %and3128, label %while.body33.preheader, label %while.cond38.preheader

while.body33.preheader:                           ; preds = %while.cond24.preheader
  br label %while.body33

while.body:                                       ; preds = %while.body.preheader, %while.body
  %p.addr.016 = phi i8* [ %add.ptr23, %while.body ], [ %p, %while.body.preheader ]
  %call5 = tail call fastcc i32 @HUF_decodeSymbolX4(i8* %p.addr.016, %struct.BIT_DStream_t* %bitDPtr, %struct.HUF_DEltX4* %dt, i32 %dtLog)
  %idx.ext = zext i32 %call5 to i64
  %add.ptr6 = getelementptr inbounds i8, i8* %p.addr.016, i64 %idx.ext
  %call10 = tail call fastcc i32 @HUF_decodeSymbolX4(i8* %add.ptr6, %struct.BIT_DStream_t* %bitDPtr, %struct.HUF_DEltX4* %dt, i32 %dtLog)
  %idx.ext11 = zext i32 %call10 to i64
  %add.ptr12 = getelementptr inbounds i8, i8* %add.ptr6, i64 %idx.ext11
  %call17 = tail call fastcc i32 @HUF_decodeSymbolX4(i8* %add.ptr12, %struct.BIT_DStream_t* %bitDPtr, %struct.HUF_DEltX4* %dt, i32 %dtLog)
  %idx.ext18 = zext i32 %call17 to i64
  %add.ptr19 = getelementptr inbounds i8, i8* %add.ptr12, i64 %idx.ext18
  %call21 = tail call fastcc i32 @HUF_decodeSymbolX4(i8* %add.ptr19, %struct.BIT_DStream_t* %bitDPtr, %struct.HUF_DEltX4* %dt, i32 %dtLog)
  %idx.ext22 = zext i32 %call21 to i64
  %add.ptr23 = getelementptr inbounds i8, i8* %add.ptr19, i64 %idx.ext22
  %call = tail call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* %bitDPtr)
  %cmp = icmp eq i32 %call, 0
  %cmp1 = icmp ult i8* %add.ptr23, %add.ptr
  %and1 = and i1 %cmp1, %cmp
  br i1 %and1, label %while.body, label %while.cond24.preheader.loopexit

while.cond38.preheader.loopexit:                  ; preds = %while.body33
  %add.ptr36.lcssa = phi i8* [ %add.ptr36, %while.body33 ]
  br label %while.cond38.preheader

while.cond38.preheader:                           ; preds = %while.cond38.preheader.loopexit, %while.cond24.preheader
  %p.addr.3.lcssa = phi i8* [ %p.addr.0.lcssa, %while.cond24.preheader ], [ %add.ptr36.lcssa, %while.cond38.preheader.loopexit ]
  %cmp403 = icmp ugt i8* %p.addr.3.lcssa, %add.ptr28
  br i1 %cmp403, label %while.end46, label %while.body42.preheader

while.body42.preheader:                           ; preds = %while.cond38.preheader
  br label %while.body42

while.body33:                                     ; preds = %while.body33.preheader, %while.body33
  %p.addr.39 = phi i8* [ %add.ptr36, %while.body33 ], [ %p.addr.0.lcssa, %while.body33.preheader ]
  %call34 = tail call fastcc i32 @HUF_decodeSymbolX4(i8* %p.addr.39, %struct.BIT_DStream_t* %bitDPtr, %struct.HUF_DEltX4* %dt, i32 %dtLog)
  %idx.ext35 = zext i32 %call34 to i64
  %add.ptr36 = getelementptr inbounds i8, i8* %p.addr.39, i64 %idx.ext35
  %call25 = tail call fastcc i32 @BIT_reloadDStream(%struct.BIT_DStream_t* %bitDPtr)
  %cmp26 = icmp eq i32 %call25, 0
  %cmp29 = icmp ule i8* %add.ptr36, %add.ptr28
  %and312 = and i1 %cmp29, %cmp26
  br i1 %and312, label %while.body33, label %while.cond38.preheader.loopexit

while.body42:                                     ; preds = %while.body42.preheader, %while.body42
  %p.addr.44 = phi i8* [ %add.ptr45, %while.body42 ], [ %p.addr.3.lcssa, %while.body42.preheader ]
  %call43 = tail call fastcc i32 @HUF_decodeSymbolX4(i8* %p.addr.44, %struct.BIT_DStream_t* %bitDPtr, %struct.HUF_DEltX4* %dt, i32 %dtLog)
  %idx.ext44 = zext i32 %call43 to i64
  %add.ptr45 = getelementptr inbounds i8, i8* %p.addr.44, i64 %idx.ext44
  %cmp40 = icmp ugt i8* %add.ptr45, %add.ptr28
  br i1 %cmp40, label %while.end46.loopexit, label %while.body42

while.end46.loopexit:                             ; preds = %while.body42
  %add.ptr45.lcssa = phi i8* [ %add.ptr45, %while.body42 ]
  br label %while.end46

while.end46:                                      ; preds = %while.end46.loopexit, %while.cond38.preheader
  %p.addr.4.lcssa = phi i8* [ %p.addr.3.lcssa, %while.cond38.preheader ], [ %add.ptr45.lcssa, %while.end46.loopexit ]
  %cmp47 = icmp ult i8* %p.addr.4.lcssa, %pEnd
  br i1 %cmp47, label %if.then49, label %if.end53

if.then49:                                        ; preds = %while.end46
  %DStream.idx.i = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitDPtr, i64 0, i32 0
  %DStream.idx.val.i = load i64, i64* %DStream.idx.i, align 8, !tbaa !40
  %DStream.idx1.i = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %bitDPtr, i64 0, i32 1
  %DStream.idx1.val.i = load i32, i32* %DStream.idx1.i, align 8, !tbaa !13
  %and.i.i = and i32 %DStream.idx1.val.i, 63
  %sh_prom.i.i = zext i32 %and.i.i to i64
  %shl.i.i = shl i64 %DStream.idx.val.i, %sh_prom.i.i
  %0 = sub i32 0, %dtLog
  %and1.i.i = and i32 %0, 63
  %sh_prom2.i.i = zext i32 %and1.i.i to i64
  %shr.i.i = lshr i64 %shl.i.i, %sh_prom2.i.i
  %add.ptr.i = getelementptr inbounds %struct.HUF_DEltX4, %struct.HUF_DEltX4* %dt, i64 %shr.i.i
  %1 = bitcast %struct.HUF_DEltX4* %add.ptr.i to i8*
  %2 = tail call i64 @llvm.objectsize.i64.p0i8(i8* %p.addr.4.lcssa, i1 false) #11
  %call1.i = tail call i8* @__memcpy_chk(i8* %p.addr.4.lcssa, i8* %1, i64 1, i64 %2) #11
  %length.i = getelementptr inbounds %struct.HUF_DEltX4, %struct.HUF_DEltX4* %dt, i64 %shr.i.i, i32 2
  %3 = load i8, i8* %length.i, align 1, !tbaa !41
  %cmp.i = icmp eq i8 %3, 1
  br i1 %cmp.i, label %if.then.i, label %if.else.i

if.then.i:                                        ; preds = %if.then49
  %nbBits.i = getelementptr inbounds %struct.HUF_DEltX4, %struct.HUF_DEltX4* %dt, i64 %shr.i.i, i32 1
  %4 = load i8, i8* %nbBits.i, align 2, !tbaa !44
  %conv4.i = zext i8 %4 to i32
  %5 = load i32, i32* %DStream.idx1.i, align 8, !tbaa !13
  %add.i3.i = add i32 %5, %conv4.i
  store i32 %add.i3.i, i32* %DStream.idx1.i, align 8, !tbaa !13
  br label %if.end53

if.else.i:                                        ; preds = %if.then49
  %6 = load i32, i32* %DStream.idx1.i, align 8, !tbaa !13
  %cmp6.i = icmp ult i32 %6, 64
  br i1 %cmp6.i, label %if.then8.i, label %if.end53

if.then8.i:                                       ; preds = %if.else.i
  %nbBits10.i = getelementptr inbounds %struct.HUF_DEltX4, %struct.HUF_DEltX4* %dt, i64 %shr.i.i, i32 1
  %7 = load i8, i8* %nbBits10.i, align 2, !tbaa !44
  %conv11.i = zext i8 %7 to i32
  %add.i.i = add i32 %conv11.i, %6
  %cmp14.i = icmp ugt i32 %add.i.i, 64
  %.add.i.i = select i1 %cmp14.i, i32 64, i32 %add.i.i
  store i32 %.add.i.i, i32* %DStream.idx1.i, align 8
  br label %if.end53

if.end53:                                         ; preds = %if.then8.i, %if.else.i, %if.then.i, %while.end46
  ret void
}

; Function Attrs: noinline nounwind ssp uwtable
define internal fastcc i32 @HUF_decodeSymbolX4(i8* %op, %struct.BIT_DStream_t* nocapture %DStream, %struct.HUF_DEltX4* %dt, i32 %dtLog) unnamed_addr #10 {
entry:
  %DStream.idx = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %DStream, i64 0, i32 0
  %DStream.idx.val = load i64, i64* %DStream.idx, align 8, !tbaa !40
  %DStream.idx11 = getelementptr %struct.BIT_DStream_t, %struct.BIT_DStream_t* %DStream, i64 0, i32 1
  %DStream.idx11.val = load i32, i32* %DStream.idx11, align 8, !tbaa !13
  %and.i = and i32 %DStream.idx11.val, 63
  %sh_prom.i = zext i32 %and.i to i64
  %shl.i = shl i64 %DStream.idx.val, %sh_prom.i
  %0 = sub i32 0, %dtLog
  %and1.i = and i32 %0, 63
  %sh_prom2.i = zext i32 %and1.i to i64
  %shr.i = lshr i64 %shl.i, %sh_prom2.i
  %add.ptr = getelementptr inbounds %struct.HUF_DEltX4, %struct.HUF_DEltX4* %dt, i64 %shr.i
  %1 = bitcast %struct.HUF_DEltX4* %add.ptr to i8*
  %2 = tail call i64 @llvm.objectsize.i64.p0i8(i8* %op, i1 false)
  %call1 = tail call i8* @__memcpy_chk(i8* %op, i8* %1, i64 2, i64 %2) #11
  %nbBits = getelementptr inbounds %struct.HUF_DEltX4, %struct.HUF_DEltX4* %dt, i64 %shr.i, i32 1
  %3 = load i8, i8* %nbBits, align 2, !tbaa !44
  %conv = zext i8 %3 to i32
  %4 = load i32, i32* %DStream.idx11, align 8, !tbaa !13
  %add.i = add i32 %4, %conv
  store i32 %add.i, i32* %DStream.idx11, align 8, !tbaa !13
  %length = getelementptr inbounds %struct.HUF_DEltX4, %struct.HUF_DEltX4* %dt, i64 %shr.i, i32 2
  %5 = load i8, i8* %length, align 1, !tbaa !41
  %conv3 = zext i8 %5 to i32
  ret i32 %conv3
}

attributes #0 = { nounwind ssp uwtable "disable-tail-calls"="false" "less-precise-fpmad"="false" "no-frame-pointer-elim"="true" "no-frame-pointer-elim-non-leaf" "no-infs-fp-math"="false" "no-nans-fp-math"="false" "stack-protector-buffer-size"="8" "target-cpu"="core2" "target-features"="+cx16,+fxsr,+mmx,+sse,+sse2,+sse3,+ssse3" "unsafe-fp-math"="false" "use-soft-float"="false" }
attributes #1 = { argmemonly nounwind }
attributes #2 = { "disable-tail-calls"="false" "less-precise-fpmad"="false" "no-frame-pointer-elim"="true" "no-frame-pointer-elim-non-leaf" "no-infs-fp-math"="false" "no-nans-fp-math"="false" "stack-protector-buffer-size"="8" "target-cpu"="core2" "target-features"="+cx16,+fxsr,+mmx,+sse,+sse2,+sse3,+ssse3" "unsafe-fp-math"="false" "use-soft-float"="false" }
attributes #3 = { nounwind "disable-tail-calls"="false" "less-precise-fpmad"="false" "no-frame-pointer-elim"="true" "no-frame-pointer-elim-non-leaf" "no-infs-fp-math"="false" "no-nans-fp-math"="false" "stack-protector-buffer-size"="8" "target-cpu"="core2" "target-features"="+cx16,+fxsr,+mmx,+sse,+sse2,+sse3,+ssse3" "unsafe-fp-math"="false" "use-soft-float"="false" }
attributes #4 = { nounwind readnone }
attributes #5 = { norecurse nounwind readnone ssp uwtable "disable-tail-calls"="false" "less-precise-fpmad"="false" "no-frame-pointer-elim"="true" "no-frame-pointer-elim-non-leaf" "no-infs-fp-math"="false" "no-nans-fp-math"="false" "stack-protector-buffer-size"="8" "target-cpu"="core2" "target-features"="+cx16,+fxsr,+mmx,+sse,+sse2,+sse3,+ssse3" "unsafe-fp-math"="false" "use-soft-float"="false" }
attributes #6 = { inlinehint nounwind ssp uwtable "disable-tail-calls"="false" "less-precise-fpmad"="false" "no-frame-pointer-elim"="true" "no-frame-pointer-elim-non-leaf" "no-infs-fp-math"="false" "no-nans-fp-math"="false" "stack-protector-buffer-size"="8" "target-cpu"="core2" "target-features"="+cx16,+fxsr,+mmx,+sse,+sse2,+sse3,+ssse3" "unsafe-fp-math"="false" "use-soft-float"="false" }
attributes #7 = { inlinehint norecurse nounwind ssp uwtable "disable-tail-calls"="false" "less-precise-fpmad"="false" "no-frame-pointer-elim"="true" "no-frame-pointer-elim-non-leaf" "no-infs-fp-math"="false" "no-nans-fp-math"="false" "stack-protector-buffer-size"="8" "target-cpu"="core2" "target-features"="+cx16,+fxsr,+mmx,+sse,+sse2,+sse3,+ssse3" "unsafe-fp-math"="false" "use-soft-float"="false" }
attributes #8 = { inlinehint noinline norecurse nounwind ssp uwtable "disable-tail-calls"="false" "less-precise-fpmad"="false" "no-frame-pointer-elim"="true" "no-frame-pointer-elim-non-leaf" "no-infs-fp-math"="false" "no-nans-fp-math"="false" "stack-protector-buffer-size"="8" "target-cpu"="core2" "target-features"="+cx16,+fxsr,+mmx,+sse,+sse2,+sse3,+ssse3" "unsafe-fp-math"="false" "use-soft-float"="false" }
attributes #9 = { noinline norecurse nounwind ssp uwtable "disable-tail-calls"="false" "less-precise-fpmad"="false" "no-frame-pointer-elim"="true" "no-frame-pointer-elim-non-leaf" "no-infs-fp-math"="false" "no-nans-fp-math"="false" "stack-protector-buffer-size"="8" "target-cpu"="core2" "target-features"="+cx16,+fxsr,+mmx,+sse,+sse2,+sse3,+ssse3" "unsafe-fp-math"="false" "use-soft-float"="false" }
attributes #10 = { noinline nounwind ssp uwtable "disable-tail-calls"="false" "less-precise-fpmad"="false" "no-frame-pointer-elim"="true" "no-frame-pointer-elim-non-leaf" "no-infs-fp-math"="false" "no-nans-fp-math"="false" "stack-protector-buffer-size"="8" "target-cpu"="core2" "target-features"="+cx16,+fxsr,+mmx,+sse,+sse2,+sse3,+ssse3" "unsafe-fp-math"="false" "use-soft-float"="false" }
attributes #11 = { nounwind }

!llvm.module.flags = !{!0}
!llvm.ident = !{!1}

!0 = !{i32 1, !"PIC Level", i32 2}
!1 = !{!"clang version 3.8.1 (tags/RELEASE_381/final)"}
!2 = !{!3, !3, i64 0}
!3 = !{!"int", !4, i64 0}
!4 = !{!"omnipotent char", !5, i64 0}
!5 = !{!"Simple C/C++ TBAA"}
!6 = !{!7, !4, i64 1}
!7 = !{!"", !4, i64 0, !4, i64 1, !4, i64 2, !4, i64 3}
!8 = !{!7, !4, i64 2}
!9 = distinct !{!9, !10}
!10 = !{!"llvm.loop.unroll.disable"}
!11 = !{!4, !4, i64 0}
!12 = distinct !{!12, !10}
!13 = !{!14, !3, i64 8}
!14 = !{!"", !15, i64 0, !3, i64 8, !16, i64 16, !16, i64 24}
!15 = !{!"long", !4, i64 0}
!16 = !{!"any pointer", !4, i64 0}
!17 = !{!14, !16, i64 16}
!18 = !{!14, !16, i64 24}
!19 = distinct !{!19, !10}
!20 = !{!21, !4, i64 0}
!21 = !{!"", !4, i64 0, !4, i64 1}
!22 = !{!21, !4, i64 1}
!23 = distinct !{!23, !10}
!24 = distinct !{!24, !25, !26}
!25 = !{!"llvm.loop.vectorize.width", i32 1}
!26 = !{!"llvm.loop.interleave.count", i32 1}
!27 = distinct !{!27, !10}
!28 = distinct !{!28, !25, !26}
!29 = distinct !{!29, !10}
!30 = distinct !{!30, !25, !26}
!31 = distinct !{!31, !32, !25, !26}
!32 = !{!"llvm.loop.unroll.runtime.disable"}
!33 = distinct !{!33, !10}
!34 = distinct !{!34, !25, !26}
!35 = distinct !{!35, !32, !25, !26}
!36 = !{!37, !3, i64 0}
!37 = !{!"", !3, i64 0, !3, i64 4}
!38 = !{!37, !3, i64 4}
!39 = !{!16, !16, i64 0}
!40 = !{!14, !15, i64 0}
!41 = !{!42, !4, i64 3}
!42 = !{!"", !43, i64 0, !4, i64 2, !4, i64 3}
!43 = !{!"short", !4, i64 0}
!44 = !{!42, !4, i64 2}
