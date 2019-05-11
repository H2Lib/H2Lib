#!/bin/sh
error=0

Tests/test_amatrix
error=`expr $error + $?` 
Tests/test_eigen
error=`expr $error + $?` 
Tests/test_hmatrix
error=`expr $error + $?` 
Tests/test_h2matrix
error=`expr $error + $?` 
Tests/test_laplacebem2d
error=`expr $error + $?` 
Tests/test_laplacebem3d
error=`expr $error + $?` 
Tests/test_laplacebem3d_ocl
error=`expr $error + $?` 
Tests/test_helmholtzbem3d
error=`expr $error + $?` 
Tests/test_helmholtzbem3d_ocl
error=`expr $error + $?` 
Tests/test_h2compression
error=`expr $error + $?` 
Tests/test_tri2dp1
error=`expr $error + $?` 
Tests/test_tet3dp1
error=`expr $error + $?` 
Tests/test_krylov
error=`expr $error + $?`
Tests/test_krylovsolvers
error=`expr $error + $?`
Tests/test_ddcluster
error=`expr $error + $?`
Tests/test_tri2drt0
error=`expr $error + $?`
Tests/test_tet3drt0
error=`expr $error + $?`
printf "\n\n  Total problems found:\n   %d\n" ${error}
