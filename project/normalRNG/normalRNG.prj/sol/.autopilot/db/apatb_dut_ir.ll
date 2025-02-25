; ModuleID = 'D:/zu7/normalRNG/normalRNG.prj/sol/.autopilot/db/a.g.ld.5.gdce.bc'
source_filename = "llvm-link"
target datalayout = "e-m:e-i64:64-i128:128-i256:256-i512:512-i1024:1024-i2048:2048-i4096:4096-n8:16:32:64-S128-v16:16-v24:32-v32:32-v48:64-v96:128-v192:256-v256:256-v512:512-v1024:1024"
target triple = "fpga64-xilinx-none"

; Function Attrs: noinline
define void @apatb_dut_ir(double* noalias nocapture nonnull readonly "fpga.decayed.dim.hint"="10" "maxi" %temp_inv, double* noalias nocapture nonnull readonly "fpga.decayed.dim.hint"="10" "maxi" %sigma, [3 x double]* noalias nocapture nonnull "fpga.decayed.dim.hint"="20000" "maxi" %sample_output, i32 %nSamples, [20000 x double]* noalias nocapture nonnull readnone "fpga.decayed.dim.hint"="10" %y, [2 x double]* noalias nocapture nonnull readonly "fpga.decayed.dim.hint"="200" %data) local_unnamed_addr #0 {
entry:
  %temp_inv_copy = alloca [10 x double], align 512
  %sigma_copy = alloca [10 x double], align 512
  %malloccall = tail call i8* @malloc(i64 480000)
  %sample_output_copy = bitcast i8* %malloccall to [20000 x [3 x double]]*
  %malloccall1 = tail call i8* @malloc(i64 1600000)
  %y_copy = bitcast i8* %malloccall1 to [10 x [20000 x double]]*
  %data_copy = alloca [200 x [2 x double]], align 512
  %0 = bitcast double* %temp_inv to [10 x double]*
  %1 = bitcast double* %sigma to [10 x double]*
  %2 = bitcast [3 x double]* %sample_output to [20000 x [3 x double]]*
  %3 = bitcast [20000 x double]* %y to [10 x [20000 x double]]*
  %4 = bitcast [2 x double]* %data to [200 x [2 x double]]*
  call fastcc void @copy_in([10 x double]* nonnull %0, [10 x double]* nonnull align 512 %temp_inv_copy, [10 x double]* nonnull %1, [10 x double]* nonnull align 512 %sigma_copy, [20000 x [3 x double]]* nonnull %2, [20000 x [3 x double]]* %sample_output_copy, [10 x [20000 x double]]* nonnull %3, [10 x [20000 x double]]* %y_copy, [200 x [2 x double]]* nonnull %4, [200 x [2 x double]]* nonnull align 512 %data_copy)
  call void @apatb_dut_hw([10 x double]* %temp_inv_copy, [10 x double]* %sigma_copy, [20000 x [3 x double]]* %sample_output_copy, i32 %nSamples, [10 x [20000 x double]]* %y_copy, [200 x [2 x double]]* %data_copy)
  call void @copy_back([10 x double]* %0, [10 x double]* %temp_inv_copy, [10 x double]* %1, [10 x double]* %sigma_copy, [20000 x [3 x double]]* %2, [20000 x [3 x double]]* %sample_output_copy, [10 x [20000 x double]]* %3, [10 x [20000 x double]]* %y_copy, [200 x [2 x double]]* %4, [200 x [2 x double]]* %data_copy)
  tail call void @free(i8* %malloccall)
  tail call void @free(i8* %malloccall1)
  ret void
}

declare noalias i8* @malloc(i64) local_unnamed_addr

; Function Attrs: argmemonly noinline norecurse willreturn
define internal fastcc void @copy_in([10 x double]* noalias readonly, [10 x double]* noalias align 512, [10 x double]* noalias readonly, [10 x double]* noalias align 512, [20000 x [3 x double]]* noalias readonly, [20000 x [3 x double]]* noalias, [10 x [20000 x double]]* noalias readonly, [10 x [20000 x double]]* noalias, [200 x [2 x double]]* noalias readonly, [200 x [2 x double]]* noalias align 512) unnamed_addr #1 {
entry:
  call fastcc void @onebyonecpy_hls.p0a10f64([10 x double]* align 512 %1, [10 x double]* %0)
  call fastcc void @onebyonecpy_hls.p0a10f64([10 x double]* align 512 %3, [10 x double]* %2)
  call fastcc void @onebyonecpy_hls.p0a20000a3f64([20000 x [3 x double]]* %5, [20000 x [3 x double]]* %4)
  call fastcc void @onebyonecpy_hls.p0a10a20000f64([10 x [20000 x double]]* %7, [10 x [20000 x double]]* %6)
  call fastcc void @onebyonecpy_hls.p0a200a2f64([200 x [2 x double]]* align 512 %9, [200 x [2 x double]]* %8)
  ret void
}

; Function Attrs: argmemonly noinline norecurse willreturn
define internal fastcc void @onebyonecpy_hls.p0a10f64([10 x double]* noalias align 512 %dst, [10 x double]* noalias readonly %src) unnamed_addr #2 {
entry:
  %0 = icmp eq [10 x double]* %dst, null
  %1 = icmp eq [10 x double]* %src, null
  %2 = or i1 %0, %1
  br i1 %2, label %ret, label %copy

copy:                                             ; preds = %entry
  call void @arraycpy_hls.p0a10f64([10 x double]* nonnull %dst, [10 x double]* nonnull %src, i64 10)
  br label %ret

ret:                                              ; preds = %copy, %entry
  ret void
}

; Function Attrs: argmemonly noinline norecurse willreturn
define void @arraycpy_hls.p0a10f64([10 x double]* %dst, [10 x double]* readonly %src, i64 %num) local_unnamed_addr #3 {
entry:
  %0 = icmp eq [10 x double]* %src, null
  %1 = icmp eq [10 x double]* %dst, null
  %2 = or i1 %1, %0
  br i1 %2, label %ret, label %copy

copy:                                             ; preds = %entry
  %for.loop.cond1 = icmp sgt i64 %num, 0
  br i1 %for.loop.cond1, label %for.loop.lr.ph, label %copy.split

for.loop.lr.ph:                                   ; preds = %copy
  br label %for.loop

for.loop:                                         ; preds = %for.loop, %for.loop.lr.ph
  %for.loop.idx2 = phi i64 [ 0, %for.loop.lr.ph ], [ %for.loop.idx.next, %for.loop ]
  %dst.addr = getelementptr [10 x double], [10 x double]* %dst, i64 0, i64 %for.loop.idx2
  %src.addr = getelementptr [10 x double], [10 x double]* %src, i64 0, i64 %for.loop.idx2
  %3 = load double, double* %src.addr, align 8
  store double %3, double* %dst.addr, align 8
  %for.loop.idx.next = add nuw nsw i64 %for.loop.idx2, 1
  %exitcond = icmp ne i64 %for.loop.idx.next, %num
  br i1 %exitcond, label %for.loop, label %copy.split

copy.split:                                       ; preds = %for.loop, %copy
  br label %ret

ret:                                              ; preds = %copy.split, %entry
  ret void
}

; Function Attrs: argmemonly noinline norecurse willreturn
define internal fastcc void @onebyonecpy_hls.p0a20000a3f64([20000 x [3 x double]]* noalias %dst, [20000 x [3 x double]]* noalias readonly %src) unnamed_addr #2 {
entry:
  %0 = icmp eq [20000 x [3 x double]]* %dst, null
  %1 = icmp eq [20000 x [3 x double]]* %src, null
  %2 = or i1 %0, %1
  br i1 %2, label %ret, label %copy

copy:                                             ; preds = %entry
  call void @arraycpy_hls.p0a20000a3f64([20000 x [3 x double]]* nonnull %dst, [20000 x [3 x double]]* nonnull %src, i64 20000)
  br label %ret

ret:                                              ; preds = %copy, %entry
  ret void
}

; Function Attrs: argmemonly noinline norecurse willreturn
define void @arraycpy_hls.p0a20000a3f64([20000 x [3 x double]]* %dst, [20000 x [3 x double]]* readonly %src, i64 %num) local_unnamed_addr #3 {
entry:
  %0 = icmp eq [20000 x [3 x double]]* %src, null
  %1 = icmp eq [20000 x [3 x double]]* %dst, null
  %2 = or i1 %1, %0
  br i1 %2, label %ret, label %copy

copy:                                             ; preds = %entry
  %for.loop.cond1 = icmp sgt i64 %num, 0
  br i1 %for.loop.cond1, label %for.loop.lr.ph, label %copy.split

for.loop.lr.ph:                                   ; preds = %copy
  br label %for.loop

for.loop:                                         ; preds = %for.loop, %for.loop.lr.ph
  %for.loop.idx2 = phi i64 [ 0, %for.loop.lr.ph ], [ %for.loop.idx.next, %for.loop ]
  %dst.addr = getelementptr [20000 x [3 x double]], [20000 x [3 x double]]* %dst, i64 0, i64 %for.loop.idx2
  %src.addr = getelementptr [20000 x [3 x double]], [20000 x [3 x double]]* %src, i64 0, i64 %for.loop.idx2
  call void @arraycpy_hls.p0a3f64([3 x double]* %dst.addr, [3 x double]* %src.addr, i64 3)
  %for.loop.idx.next = add nuw nsw i64 %for.loop.idx2, 1
  %exitcond = icmp ne i64 %for.loop.idx.next, %num
  br i1 %exitcond, label %for.loop, label %copy.split

copy.split:                                       ; preds = %for.loop, %copy
  br label %ret

ret:                                              ; preds = %copy.split, %entry
  ret void
}

; Function Attrs: argmemonly noinline norecurse willreturn
define void @arraycpy_hls.p0a3f64([3 x double]* %dst, [3 x double]* readonly %src, i64 %num) local_unnamed_addr #3 {
entry:
  %0 = icmp eq [3 x double]* %src, null
  %1 = icmp eq [3 x double]* %dst, null
  %2 = or i1 %1, %0
  br i1 %2, label %ret, label %copy

copy:                                             ; preds = %entry
  %for.loop.cond1 = icmp sgt i64 %num, 0
  br i1 %for.loop.cond1, label %for.loop.lr.ph, label %copy.split

for.loop.lr.ph:                                   ; preds = %copy
  br label %for.loop

for.loop:                                         ; preds = %for.loop, %for.loop.lr.ph
  %for.loop.idx2 = phi i64 [ 0, %for.loop.lr.ph ], [ %for.loop.idx.next, %for.loop ]
  %dst.addr = getelementptr [3 x double], [3 x double]* %dst, i64 0, i64 %for.loop.idx2
  %src.addr = getelementptr [3 x double], [3 x double]* %src, i64 0, i64 %for.loop.idx2
  %3 = load double, double* %src.addr, align 8
  store double %3, double* %dst.addr, align 8
  %for.loop.idx.next = add nuw nsw i64 %for.loop.idx2, 1
  %exitcond = icmp ne i64 %for.loop.idx.next, %num
  br i1 %exitcond, label %for.loop, label %copy.split

copy.split:                                       ; preds = %for.loop, %copy
  br label %ret

ret:                                              ; preds = %copy.split, %entry
  ret void
}

; Function Attrs: argmemonly noinline norecurse willreturn
define internal fastcc void @onebyonecpy_hls.p0a10a20000f64([10 x [20000 x double]]* noalias %dst, [10 x [20000 x double]]* noalias readonly %src) unnamed_addr #2 {
entry:
  %0 = icmp eq [10 x [20000 x double]]* %dst, null
  %1 = icmp eq [10 x [20000 x double]]* %src, null
  %2 = or i1 %0, %1
  br i1 %2, label %ret, label %copy

copy:                                             ; preds = %entry
  call void @arraycpy_hls.p0a10a20000f64([10 x [20000 x double]]* nonnull %dst, [10 x [20000 x double]]* nonnull %src, i64 10)
  br label %ret

ret:                                              ; preds = %copy, %entry
  ret void
}

; Function Attrs: argmemonly noinline norecurse willreturn
define void @arraycpy_hls.p0a10a20000f64([10 x [20000 x double]]* %dst, [10 x [20000 x double]]* readonly %src, i64 %num) local_unnamed_addr #3 {
entry:
  %0 = icmp eq [10 x [20000 x double]]* %src, null
  %1 = icmp eq [10 x [20000 x double]]* %dst, null
  %2 = or i1 %1, %0
  br i1 %2, label %ret, label %copy

copy:                                             ; preds = %entry
  %for.loop.cond1 = icmp sgt i64 %num, 0
  br i1 %for.loop.cond1, label %for.loop.lr.ph, label %copy.split

for.loop.lr.ph:                                   ; preds = %copy
  br label %for.loop

for.loop:                                         ; preds = %for.loop, %for.loop.lr.ph
  %for.loop.idx2 = phi i64 [ 0, %for.loop.lr.ph ], [ %for.loop.idx.next, %for.loop ]
  %dst.addr = getelementptr [10 x [20000 x double]], [10 x [20000 x double]]* %dst, i64 0, i64 %for.loop.idx2
  %src.addr = getelementptr [10 x [20000 x double]], [10 x [20000 x double]]* %src, i64 0, i64 %for.loop.idx2
  call void @arraycpy_hls.p0a20000f64([20000 x double]* %dst.addr, [20000 x double]* %src.addr, i64 20000)
  %for.loop.idx.next = add nuw nsw i64 %for.loop.idx2, 1
  %exitcond = icmp ne i64 %for.loop.idx.next, %num
  br i1 %exitcond, label %for.loop, label %copy.split

copy.split:                                       ; preds = %for.loop, %copy
  br label %ret

ret:                                              ; preds = %copy.split, %entry
  ret void
}

; Function Attrs: argmemonly noinline norecurse willreturn
define void @arraycpy_hls.p0a20000f64([20000 x double]* %dst, [20000 x double]* readonly %src, i64 %num) local_unnamed_addr #3 {
entry:
  %0 = icmp eq [20000 x double]* %src, null
  %1 = icmp eq [20000 x double]* %dst, null
  %2 = or i1 %1, %0
  br i1 %2, label %ret, label %copy

copy:                                             ; preds = %entry
  %for.loop.cond1 = icmp sgt i64 %num, 0
  br i1 %for.loop.cond1, label %for.loop.lr.ph, label %copy.split

for.loop.lr.ph:                                   ; preds = %copy
  br label %for.loop

for.loop:                                         ; preds = %for.loop, %for.loop.lr.ph
  %for.loop.idx2 = phi i64 [ 0, %for.loop.lr.ph ], [ %for.loop.idx.next, %for.loop ]
  %dst.addr = getelementptr [20000 x double], [20000 x double]* %dst, i64 0, i64 %for.loop.idx2
  %src.addr = getelementptr [20000 x double], [20000 x double]* %src, i64 0, i64 %for.loop.idx2
  %3 = load double, double* %src.addr, align 8
  store double %3, double* %dst.addr, align 8
  %for.loop.idx.next = add nuw nsw i64 %for.loop.idx2, 1
  %exitcond = icmp ne i64 %for.loop.idx.next, %num
  br i1 %exitcond, label %for.loop, label %copy.split

copy.split:                                       ; preds = %for.loop, %copy
  br label %ret

ret:                                              ; preds = %copy.split, %entry
  ret void
}

; Function Attrs: argmemonly noinline norecurse willreturn
define internal fastcc void @onebyonecpy_hls.p0a200a2f64([200 x [2 x double]]* noalias align 512 %dst, [200 x [2 x double]]* noalias readonly %src) unnamed_addr #2 {
entry:
  %0 = icmp eq [200 x [2 x double]]* %dst, null
  %1 = icmp eq [200 x [2 x double]]* %src, null
  %2 = or i1 %0, %1
  br i1 %2, label %ret, label %copy

copy:                                             ; preds = %entry
  call void @arraycpy_hls.p0a200a2f64([200 x [2 x double]]* nonnull %dst, [200 x [2 x double]]* nonnull %src, i64 200)
  br label %ret

ret:                                              ; preds = %copy, %entry
  ret void
}

; Function Attrs: argmemonly noinline norecurse willreturn
define void @arraycpy_hls.p0a200a2f64([200 x [2 x double]]* %dst, [200 x [2 x double]]* readonly %src, i64 %num) local_unnamed_addr #3 {
entry:
  %0 = icmp eq [200 x [2 x double]]* %src, null
  %1 = icmp eq [200 x [2 x double]]* %dst, null
  %2 = or i1 %1, %0
  br i1 %2, label %ret, label %copy

copy:                                             ; preds = %entry
  %for.loop.cond1 = icmp sgt i64 %num, 0
  br i1 %for.loop.cond1, label %for.loop.lr.ph, label %copy.split

for.loop.lr.ph:                                   ; preds = %copy
  br label %for.loop

for.loop:                                         ; preds = %for.loop, %for.loop.lr.ph
  %for.loop.idx2 = phi i64 [ 0, %for.loop.lr.ph ], [ %for.loop.idx.next, %for.loop ]
  %dst.addr = getelementptr [200 x [2 x double]], [200 x [2 x double]]* %dst, i64 0, i64 %for.loop.idx2
  %src.addr = getelementptr [200 x [2 x double]], [200 x [2 x double]]* %src, i64 0, i64 %for.loop.idx2
  call void @arraycpy_hls.p0a2f64([2 x double]* %dst.addr, [2 x double]* %src.addr, i64 2)
  %for.loop.idx.next = add nuw nsw i64 %for.loop.idx2, 1
  %exitcond = icmp ne i64 %for.loop.idx.next, %num
  br i1 %exitcond, label %for.loop, label %copy.split

copy.split:                                       ; preds = %for.loop, %copy
  br label %ret

ret:                                              ; preds = %copy.split, %entry
  ret void
}

; Function Attrs: argmemonly noinline norecurse willreturn
define void @arraycpy_hls.p0a2f64([2 x double]* %dst, [2 x double]* readonly %src, i64 %num) local_unnamed_addr #3 {
entry:
  %0 = icmp eq [2 x double]* %src, null
  %1 = icmp eq [2 x double]* %dst, null
  %2 = or i1 %1, %0
  br i1 %2, label %ret, label %copy

copy:                                             ; preds = %entry
  %for.loop.cond1 = icmp sgt i64 %num, 0
  br i1 %for.loop.cond1, label %for.loop.lr.ph, label %copy.split

for.loop.lr.ph:                                   ; preds = %copy
  br label %for.loop

for.loop:                                         ; preds = %for.loop, %for.loop.lr.ph
  %for.loop.idx2 = phi i64 [ 0, %for.loop.lr.ph ], [ %for.loop.idx.next, %for.loop ]
  %dst.addr = getelementptr [2 x double], [2 x double]* %dst, i64 0, i64 %for.loop.idx2
  %src.addr = getelementptr [2 x double], [2 x double]* %src, i64 0, i64 %for.loop.idx2
  %3 = load double, double* %src.addr, align 8
  store double %3, double* %dst.addr, align 8
  %for.loop.idx.next = add nuw nsw i64 %for.loop.idx2, 1
  %exitcond = icmp ne i64 %for.loop.idx.next, %num
  br i1 %exitcond, label %for.loop, label %copy.split

copy.split:                                       ; preds = %for.loop, %copy
  br label %ret

ret:                                              ; preds = %copy.split, %entry
  ret void
}

; Function Attrs: argmemonly noinline norecurse willreturn
define internal fastcc void @copy_out([10 x double]* noalias, [10 x double]* noalias readonly align 512, [10 x double]* noalias, [10 x double]* noalias readonly align 512, [20000 x [3 x double]]* noalias, [20000 x [3 x double]]* noalias readonly, [10 x [20000 x double]]* noalias, [10 x [20000 x double]]* noalias readonly, [200 x [2 x double]]* noalias, [200 x [2 x double]]* noalias readonly align 512) unnamed_addr #4 {
entry:
  call fastcc void @onebyonecpy_hls.p0a10f64([10 x double]* %0, [10 x double]* align 512 %1)
  call fastcc void @onebyonecpy_hls.p0a10f64([10 x double]* %2, [10 x double]* align 512 %3)
  call fastcc void @onebyonecpy_hls.p0a20000a3f64([20000 x [3 x double]]* %4, [20000 x [3 x double]]* %5)
  call fastcc void @onebyonecpy_hls.p0a10a20000f64([10 x [20000 x double]]* %6, [10 x [20000 x double]]* %7)
  call fastcc void @onebyonecpy_hls.p0a200a2f64([200 x [2 x double]]* %8, [200 x [2 x double]]* align 512 %9)
  ret void
}

declare void @free(i8*) local_unnamed_addr

declare void @apatb_dut_hw([10 x double]*, [10 x double]*, [20000 x [3 x double]]*, i32, [10 x [20000 x double]]*, [200 x [2 x double]]*)

; Function Attrs: argmemonly noinline norecurse willreturn
define internal fastcc void @copy_back([10 x double]* noalias, [10 x double]* noalias readonly align 512, [10 x double]* noalias, [10 x double]* noalias readonly align 512, [20000 x [3 x double]]* noalias, [20000 x [3 x double]]* noalias readonly, [10 x [20000 x double]]* noalias, [10 x [20000 x double]]* noalias readonly, [200 x [2 x double]]* noalias, [200 x [2 x double]]* noalias readonly align 512) unnamed_addr #4 {
entry:
  call fastcc void @onebyonecpy_hls.p0a20000a3f64([20000 x [3 x double]]* %4, [20000 x [3 x double]]* %5)
  ret void
}

define void @dut_hw_stub_wrapper([10 x double]*, [10 x double]*, [20000 x [3 x double]]*, i32, [10 x [20000 x double]]*, [200 x [2 x double]]*) #5 {
entry:
  call void @copy_out([10 x double]* null, [10 x double]* %0, [10 x double]* null, [10 x double]* %1, [20000 x [3 x double]]* null, [20000 x [3 x double]]* %2, [10 x [20000 x double]]* null, [10 x [20000 x double]]* %4, [200 x [2 x double]]* null, [200 x [2 x double]]* %5)
  %6 = bitcast [10 x double]* %0 to double*
  %7 = bitcast [10 x double]* %1 to double*
  %8 = bitcast [20000 x [3 x double]]* %2 to [3 x double]*
  %9 = bitcast [10 x [20000 x double]]* %4 to [20000 x double]*
  %10 = bitcast [200 x [2 x double]]* %5 to [2 x double]*
  call void @dut_hw_stub(double* %6, double* %7, [3 x double]* %8, i32 %3, [20000 x double]* %9, [2 x double]* %10)
  call void @copy_in([10 x double]* null, [10 x double]* %0, [10 x double]* null, [10 x double]* %1, [20000 x [3 x double]]* null, [20000 x [3 x double]]* %2, [10 x [20000 x double]]* null, [10 x [20000 x double]]* %4, [200 x [2 x double]]* null, [200 x [2 x double]]* %5)
  ret void
}

declare void @dut_hw_stub(double*, double*, [3 x double]*, i32, [20000 x double]*, [2 x double]*)

attributes #0 = { noinline "fpga.wrapper.func"="wrapper" }
attributes #1 = { argmemonly noinline norecurse willreturn "fpga.wrapper.func"="copyin" }
attributes #2 = { argmemonly noinline norecurse willreturn "fpga.wrapper.func"="onebyonecpy_hls" }
attributes #3 = { argmemonly noinline norecurse willreturn "fpga.wrapper.func"="arraycpy_hls" }
attributes #4 = { argmemonly noinline norecurse willreturn "fpga.wrapper.func"="copyout" }
attributes #5 = { "fpga.wrapper.func"="stub" }

!llvm.dbg.cu = !{}
!llvm.ident = !{!0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0, !0}
!llvm.module.flags = !{!1, !2, !3}
!blackbox_cfg = !{!4}

!0 = !{!"clang version 7.0.0 "}
!1 = !{i32 2, !"Dwarf Version", i32 4}
!2 = !{i32 2, !"Debug Info Version", i32 3}
!3 = !{i32 1, !"wchar_size", i32 4}
!4 = !{}
