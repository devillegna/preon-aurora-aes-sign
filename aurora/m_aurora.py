
import sys
sys.path.insert(0, '../' )

import time

from libgf264fft import clib_wrapper as cgf

# defining gf
from utils import gf2192 as gf

from utils import randombytes as rd

# the commit scheme
from commit import merkeltreecommit as mt

# hash for challenges
from utils import hash as H

from fri_ldt import m_fri_ldt as fri

import numpy as np

def _dummy( *params ) : return None

def _log2( w ):
    assert 0 != w
    assert 0== (w&(w-1))
    r, t = 0 , 1
    while(w!=t): r , t = r+1 , t<<1
    return r

def _pad_len( n_or_m ):
    r = 1
    while( r < n_or_m ):  r = r<<1
    return r

def commit_polys( polys , max_poly_len , RS_rho = 8 , RS_shift = (1<<63) , verbose = 1 ):
    if 1 == verbose : dump = print
    else : dump = _dummy

    dump( "commit polys: (poly_len):" , [ len(p) for p in polys ] )
    mesg0 = []
    for pp in polys :
        v_p0 = gf.fft( pp , RS_rho * (max_poly_len//len(pp)) , RS_shift )
        v_p1 = [ int.to_bytes(x,gf.GF_BSIZE,'little') for x in v_p0 ]
        v_p2 = [ b''.join( (v_p1[2*i] , v_p1[2*i+1]) ) for i in range(len(v_p1)//2) ]
        mesg0.append( v_p2 )
    dump( f"mesg0: ({max_poly_len}x{RS_rho}/2 = {max_poly_len*RS_rho//2}) " , [len(v) for v in mesg0] )
    mesgs = [ b''.join(e) for e in zip(*mesg0) ]
    dump( f"|mesgs:| = {len(mesgs)} , |mesgs[0]| = {len(mesgs[0])} == {gf.GF_BSIZE * 2 * len(mesg0)}" )
    rt , r_leaf , mktree = mt.commit( mesgs )
    return rt , mesgs , r_leaf , mktree

def mat_x_vec( mat , vec , n_row=-1) :
    if 0 > n_row : n_row = len(mat)
    r = [0]*n_row
    for i in range(n_row):
        ri = 0
        for j in range(len(vec)):
            if 0 == mat[i][j] : continue     # it's ok since the mat is public.
            ri ^= cgf.gf264_mul( mat[i][j] , vec[j] )
        r[i] = ri
    return r

def process_R1CS( R1CS , verbose = 1 ):
    if 1 == verbose : dump = print
    else : dump = _dummy

    mat_A , mat_B , mat_C , vec_z , witness_idx = R1CS
    n = len(vec_z)
    m = len(mat_A)
    pad_len = _pad_len( max(n,m) )
    inst_dim = _log2(witness_idx)
    dump( f"pad_len: {pad_len}, inst_dim: {inst_dim}" )

    dump( "padding and calcuate f_1v, f_z, f_w" )
    st = time.time()
    p_vec_z = vec_z + [ gf.from_bytes(rd.randombytes(gf.GF_BSIZE)) for _ in range(pad_len-len(vec_z)) ]   # ZK
    #p_vec_z = vec_z + [ 0 for _ in range(pad_len-len(vec_z)) ]    # no ZK
    f_z = gf.ifft( p_vec_z , 1 , 0 );   p_vec_z.extend( gf.fft( f_z , 1 , pad_len ) )
    f_w = gf.polydiv( f_z , inst_dim )[witness_idx:] + [0]*witness_idx
    ed = time.time()
    dump( "time:", format(ed-st) , "secs" )
    dump( f"f_z : [{len(f_z)}] ...[{witness_idx-4}:{witness_idx+4}]:", f_z[witness_idx-4:witness_idx+4] , "..." )
    dump( f"f_w : (f_z - f_1v) /Z_{inst_dim} : " , f_w[:4] , "..." )

    #f_1v = cgf.ibtfy( vec_z[:witness_idx] , 0 )   # shouldn't compute it actually. do this for checking correctness only.
    #dump( f"f_1v : [{len(f_1v)}] ...[{witness_idx-4}:]:", f_1v[witness_idx-4:] )

    dump( "calculate Az, Bz, Cz" )
    st = time.time()
    Az = mat_x_vec( mat_A , vec_z )
    Bz = mat_x_vec( mat_B , vec_z )
    Cz = mat_x_vec( mat_C , vec_z )
    ed = time.time()
    dump( "time:", format(ed-st) , "secs" )
    # XXX: check if Az*Bz == Cz here

    dump( "calculate f_Az, f_Bz, f_Cz" )
    st = time.time()
    v_Az    = Az + [ gf.from_bytes(rd.randombytes(gf.GF_BSIZE)) for _ in range(pad_len-len(Az)) ]
    v_Bz    = Bz + [ gf.from_bytes(rd.randombytes(gf.GF_BSIZE)) for _ in range(pad_len-len(Bz)) ]
    v_Cz    = Cz + [ gf.mul(v_Az[i],v_Bz[i]) for i in range(m,pad_len) ]           # !!! XXX: need to discuss here
    f_Az = gf.ifft( v_Az , 1 , 0 )
    f_Bz = gf.ifft( v_Bz , 1 , 0 )
    f_Cz = gf.ifft( v_Cz , 1 , 0 )
    v_Az.extend( gf.fft( f_Az , 1 , pad_len ) )
    v_Bz.extend( gf.fft( f_Bz , 1 , pad_len ) )
    v_Cz.extend( gf.fft( f_Cz , 1 , pad_len ) )
    ed = time.time()
    dump( "time:", format(ed-st) , "secs" )

    return f_w , p_vec_z , f_Az , f_Bz , f_Cz , v_Az , v_Bz , v_Cz

def row_check( v_Az , v_Bz , v_Cz , pad_len , verbose = 1 ) :
    if 1 == verbose : dump = print
    else : dump = _dummy

    ## row check
    dump( f"row check: calculate (f_Azxf_Bz - f_Cz) / Z_{pad_len}" )
    st = time.time()
    v_AzxBz = v_Cz[:pad_len] + [ gf.mul(v_Az[i],v_Bz[i]) for i in range(pad_len,2*pad_len) ]
    f_Azxf_Bz = gf.ifft( v_AzxBz , 1 , 0 )
    ed = time.time()
    dump( "time:", format(ed-st) , "secs" )
    dump( f"f_Azxf_Bz: [{len(f_Azxf_Bz)}] ...[{pad_len-2}:{pad_len+2}]" , f_Azxf_Bz[pad_len-2:pad_len+2] , "..." )
    #dump( f"f_Cz     : [{len(f_Cz)}] ...[{pad_len-4}:]" , f_Cz[pad_len-4:] )
    f_rowcheck = f_Azxf_Bz[pad_len:]
    dump( f"f_AB-C / Zh: [{len(f_rowcheck)}] " , f_rowcheck[:2] , "..." )
    return f_rowcheck

def lincheck_step1( alpha , mat_A , mat_B , mat_C , pad_len , value_or_poly , verbose = 1 ) :
    if 1 == verbose : dump = print
    else : dump = _dummy

    m = len(mat_A)
    n = len(mat_A[0])

    ## lincheck
    dump( "lin-check starts. calculate f_alpha, v_alpha" )
    st = time.time()
    v_alpha =  [ 1 , alpha ] + [0]*(pad_len-2)
    for i in range(2,m): v_alpha[i] = gf.mul( v_alpha[i-1], v_alpha[1] )
    f_alpha = gf.ifft( v_alpha , 1 , 0 )
    ed = time.time()
    dump( "time:" , format(ed-st) , "secs" )

    dump( "lin-check step 1. calculate p2A, p2B, p2C and evaluate their values" )
    dump( "evaluate values of p2A, p2B, p2C" )
    st = time.time()
    v_a_264 = list( zip( *[gf.to_gf264s(e) for e in v_alpha[:m]] ) )
    #v_a_264 = [ np.array(e,dtype=np.uint64) for e in _v_a_264 ]
    _v_p2A_264 = []
    for la in v_a_264 :
          _v_p2Ai = [0]*n
          for i in range(m):
              for j in range(n) :
                  if(0==mat_A[i][j]) : continue
                  _v_p2Ai[j] ^= cgf.gf264_mul(la[i],mat_A[i][j])
          _v_p2A_264.append( _v_p2Ai )
    v_p2A = [ gf.from_gf264s(*e) for e in zip(*_v_p2A_264) ] + [0]*(pad_len-n)
    _v_p2B_264 = []
    for la in v_a_264 :
          _v_p2Bi = [0]*n
          for i in range(m):
              for j in range(n) :
                  if(0==mat_B[i][j]) : continue
                  _v_p2Bi[j] ^= cgf.gf264_mul(la[i],mat_B[i][j])
          _v_p2B_264.append( _v_p2Bi )
    v_p2B = [ gf.from_gf264s(*e) for e in zip(*_v_p2B_264) ] + [0]*(pad_len-n)
    _v_p2C_264 = []
    for la in v_a_264 :
          _v_p2Ci = [0]*n
          for i in range(m):
              for j in range(n) :
                  if(0==mat_C[i][j]) : continue
                  _v_p2Ci[j] ^= cgf.gf264_mul(la[i],mat_C[i][j])
          _v_p2C_264.append( _v_p2Ci )
    v_p2C = [ gf.from_gf264s(*e) for e in zip(*_v_p2C_264) ] + [0]*(pad_len-n)
    ed = time.time()
    dump( "time:" , format(ed-st) , "secs" )
    # simple but slow code
    #st = time.time()
    #v_p2A , v_p2B , v_p2C = [0]*pad_len , [0]*pad_len , [0]*pad_len
    #for j in range(n):
    #    aj, bj, cj = mat_A[0][j] , mat_B[0][j] , mat_C[0][j]
    #    for i in range(1,m):
    #        aj ^= gf.mul_gf264( v_alpha[i] , mat_A[i][j] )
    #        bj ^= gf.mul_gf264( v_alpha[i] , mat_B[i][j] )
    #        cj ^= gf.mul_gf264( v_alpha[i] , mat_C[i][j] )
    #    v_p2A[j] , v_p2B[j] , v_p2C[j] = aj , bj , cj
    #ed = time.time()
    #dump( "time:" , format(ed-st) , "secs" )
    dump( "interpolate p2A, p2B, p2C" )
    st = time.time()
    p2A , p2B , p2C = gf.ifft( v_p2A , 1 , 0 ) , gf.ifft( v_p2B , 1 , 0 ) , gf.ifft( v_p2C , 1 , 0 )
    if value_or_poly :
        v_alpha.extend( gf.fft(f_alpha,1,pad_len) )
        v_p2A.extend( gf.fft(p2A,1,pad_len) )
        v_p2B.extend( gf.fft(p2B,1,pad_len) )
        v_p2C.extend( gf.fft(p2C,1,pad_len) )
    ed = time.time()
    dump( "time:" , format(ed-st) , "secs" )
    dump( "return v_alpha, v_p2A, v_p2B, b_p2C" if value_or_poly else "return f_alpha, p2A , p2B , p2C" )
    return (v_alpha , v_p2A , v_p2B , v_p2C ) if value_or_poly else ( f_alpha, p2A , p2B , p2C )


def lincheck_step2( v_alpha , v_p2A , v_p2B , v_p2C ,  p_vec_z , v_Az , v_Bz , v_Cz ,
              s1 , s2 , s3 , r_lincheck , pad_len , verbose = 1 ) :
    if 1 == verbose : dump = print
    else : dump = _dummy

    dump( f"lin-check step 2. poly muls and /Z_{_log2(pad_len)}" )
    st = time.time()
    v_sA = [ gf.mul(v_Az[i],v_alpha[i]) ^ gf.mul(p_vec_z[i],v_p2A[i]) for i in range(2*pad_len) ]
    v_sB = [ gf.mul(v_Bz[i],v_alpha[i]) ^ gf.mul(p_vec_z[i],v_p2B[i]) for i in range(2*pad_len) ]
    v_sC = [ gf.mul(v_Cz[i],v_alpha[i]) ^ gf.mul(p_vec_z[i],v_p2C[i]) for i in range(2*pad_len) ]
    f_sA = gf.ifft( v_sA , 1 , 0 )
    f_sB = gf.ifft( v_sB , 1 , 0 )
    f_sC = gf.ifft( v_sC , 1 , 0 )
    g = [ gf.mul(s1,f_sA[i])^gf.mul(s2,f_sB[i])^gf.mul(s3,f_sC[i])^r_lincheck[i] for i in range(pad_len) ]
    h = [ gf.mul(s1,f_sA[i])^gf.mul(s2,f_sB[i])^gf.mul(s3,f_sC[i])^r_lincheck[i] for i in range(pad_len,2*pad_len) ]
    ed = time.time()
    dump( "time:" , format(ed-st) , "secs" )
    dump( f"g: [{pad_len}] ...[{pad_len-2}:{pad_len}]", g[pad_len-2:] )
    dump( f"h: [{pad_len}] ...[{pad_len-2}:{pad_len}]", h[pad_len-2:] )

    return g , h

def generate_proof( R1CS , h_state , Nq = 26 , RS_rho = 8 , RS_shift=1<<63 , verbose = 1 ) :
    if 1 == verbose : dump = print
    else : dump = _dummy

    ## process R1CS
    mat_A , mat_B , mat_C , vec_z , witness_idx = R1CS
    n = len(vec_z)
    m = len(mat_A)
    pad_len = _pad_len( max(n,m) )
    dump( f"m(#rows): {m} x n: {n}, witness_idx: {witness_idx}, pad_len: {pad_len}" )

    f_w ,  p_vec_z , f_Az , f_Bz , f_Cz , v_Az , v_Bz , v_Cz = process_R1CS( R1CS , verbose )

    ## row check
    dump( "rowcheck" )
    f_rowcheck = row_check(v_Az, v_Bz, v_Cz , pad_len , verbose=0 )

    ## generate random polynomials of degree 2xpad_len: r_lincheck and r_ldt
    v_r_lincheck = [ gf.from_bytes(rd.randombytes(gf.GF_BSIZE)) for _ in range(pad_len*2) ]
    for i in range(pad_len): v_r_lincheck[0] ^= v_r_lincheck[i]
    r_lincheck = gf.ifft(v_r_lincheck, 1 , 0 )
    dump( f"r_lincheck: [{len(r_lincheck)}]: ...[{pad_len-2}:{pad_len+2}] ...", r_lincheck[pad_len-2:pad_len+2] )
    r_ldt = [ gf.from_bytes(rd.randombytes(gf.GF_BSIZE)) for _ in range(pad_len*2) ]

    proof = []
    ## first commit them: HASH their RS codeword
    dump( "commit f_w, f_Az, f_Bz , f_Cz, r_lincheck , r_ldt" )
    rt0 , mesgs0 , r_leaf0 , mktree0 = commit_polys( [f_w , f_Az , f_Bz , f_Cz , r_lincheck , r_ldt ] , 2*pad_len , RS_rho , RS_shift , verbose )
    proof.append( rt0 )
    h_state = H.gen( h_state , rt0 )

    ## lin check
    chals = [ H.gen( h_state , bytes([1,i]) )[:gf.GF_BSIZE] for i in range(1,5) ]
    alpha, s1, s2, s3 = gf.from_bytes(chals[0]), gf.from_bytes(chals[1]), gf.from_bytes(chals[2]), gf.from_bytes(chals[3])
    vs    = lincheck_step1( alpha , mat_A , mat_B , mat_C , pad_len , 1 , verbose )
    g , h = lincheck_step2( *vs , p_vec_z , v_Az , v_Bz , v_Cz , s1 , s2 , s3  , r_lincheck , pad_len , verbose )
    dump( f"commit h" )
    rt1 , mesgs1 , r_leaf1 , mktree1 = commit_polys( [ h ] , 2*pad_len , RS_rho , RS_shift , verbose )
    proof.append( rt1 )
    h_state = H.gen( *chals )

    ## generate f0 for fri_ldt and perform fri_ldt
    y = [ gf.from_bytes( H.gen( h_state , bytes([2,i]) )[:gf.GF_BSIZE] ) for i in range(1,10) ]
    g_raise = gf.ipolydiv([0]+g[:pad_len-1],0)  # raise g by degree 1. It still has to be raised by pad_len.
    f0 = [ gf.mul(y[0],f_w[i])^gf.mul(y[1],f_Az[i])^gf.mul(y[2],f_Bz[i])^gf.mul(y[3],f_Cz[i])
          ^gf.mul(y[4],f_rowcheck[i])
          ^gf.mul(y[5],r_lincheck[i])^gf.mul(y[6],h[i]) ^r_ldt[i]
          ^gf.mul(y[7],g[i])
           for i in range(pad_len) ] + [ 
           gf.mul(y[5],r_lincheck[pad_len+i])^r_ldt[pad_len+i]
          ^gf.mul(y[8],rgi)
           for i,rgi in enumerate(g_raise) ]
    ## LDT f0
    dump( "ldt |f0|:", len(f0) )
    v_f0 = gf.fft( f0 , RS_rho , RS_shift )
    dump( "calculate RS code of f0: |v_f0|: " , len(v_f0) )
    ldt_commits , ldt_d1poly , ldt_mktrees , h_state = fri.ldt_commit_phase( v_f0 , len(f0) , h_state , RS_rho , RS_shift, verbose=0 )
    ldt_open_mesgs , ldt_queries = fri.ldt_query_phase( len(f0) , ldt_mktrees , h_state , Nq , RS_rho , verbose=0 )
    dump( "ldt queries:" , ldt_queries )

    proof.extend( ldt_commits )
    proof.append( ldt_d1poly )
    proof.extend( ldt_open_mesgs )

    ## open queries
    proof.append( mt.batchopen(ldt_queries,mesgs0,r_leaf0,mktree0) )
    proof.append( mt.batchopen(ldt_queries,mesgs1,r_leaf1,mktree1) )

    return proof

def codewords_of_public_polynomials( alpha , vec_1v , mat_A , mat_B , mat_C , pad_len , RS_rho = 8 , RS_shift = (1<<63) , verbose = 1 ):
    if 1 == verbose : dump = print
    else : dump = _dummy
    dump( "prepare f_alpha, p2A, p2B, p2C" )
    f_alpha , p2A , p2B , p2C = lincheck_step1( alpha , mat_A , mat_B , mat_C , pad_len , 0 , verbose )
    dump( "generate RS codeword of public polynomials" )
    f_1v = cgf.ibtfy( vec_1v , 0 )
    rs_f_1v    = gf.fft( f_1v , (pad_len//(len(f_1v)))*2*RS_rho , RS_shift )
    rs_f_alpha = gf.fft( f_alpha , 2*RS_rho , RS_shift )
    rs_f_p2A = gf.fft( p2A , 2*RS_rho , RS_shift )
    rs_f_p2B = gf.fft( p2B , 2*RS_rho , RS_shift )
    rs_f_p2C = gf.fft( p2C , 2*RS_rho , RS_shift )
    return rs_f_1v, rs_f_alpha, rs_f_p2A, rs_f_p2B, rs_f_p2C

def values_from_virtual_oracle( _idx , aurora_open0 , aurora_open1 , lincheck_s , y , rs_codewords , inst_dim , r1cs_dim , RS_shift = (1<<63) ):
    offset = RS_shift
    # unpack input
    rs_f_1v , rs_f_alpha , rs_f_p2A, rs_f_p2B, rs_f_p2C = rs_codewords
    s1 , s2, s3 = lincheck_s
    vv0 = [ gf.from_bytes_x2(aurora_open0[i*gf.GF_BSIZE*2:i*gf.GF_BSIZE*2+gf.GF_BSIZE*2]) for i in range(6) ]
    v_w0 , v_Az0 , v_Bz0 , v_Cz0 , v_lincheck0 , v_ldt0 = vv0[0],vv0[1],vv0[2],vv0[3],vv0[4],vv0[5]
    v_h0   = gf.from_bytes_x2(aurora_open1)

    # generate output
    values = []
    for i in range(2):
        idx = _idx*2 + i
        cc0  = gf.mul(y[0],v_w0[i])^gf.mul(y[1],v_Az0[i])^gf.mul(y[2],v_Bz0[i])^gf.mul(y[3],v_Cz0[i])
        v_f_rowcheck0 = gf.mul( (gf.mul( v_Az0[i] , v_Bz0[i] )^v_Cz0[i]) , cgf.gf264_inv( cgf.index_to_gf264((offset+idx)>>r1cs_dim) ) )
        cc0 ^= gf.mul(y[4],v_f_rowcheck0)
        cc0 ^= gf.mul(y[5],v_lincheck0[i])^gf.mul(y[6],v_h0[i])^v_ldt0[i]
    
        v_fz0 = rs_f_1v[idx]^gf.mul_gf264( v_w0[i] , cgf.index_to_gf264( (offset+idx)>>inst_dim ) )
    
        v_g0 =  gf.mul( s1 , gf.mul(v_Az0[i],rs_f_alpha[idx])^gf.mul(rs_f_p2A[idx],v_fz0) ) \
               ^gf.mul( s2 , gf.mul(v_Bz0[i],rs_f_alpha[idx])^gf.mul(rs_f_p2B[idx],v_fz0) ) \
               ^gf.mul( s3 , gf.mul(v_Cz0[i],rs_f_alpha[idx])^gf.mul(rs_f_p2C[idx],v_fz0) ) \
               ^ v_lincheck0[i] ^ gf.mul_gf264( v_h0[i] , cgf.index_to_gf264((offset+idx)>>r1cs_dim) )
        cc0 ^= gf.mul(y[7],v_g0) 
        cc0 ^= gf.mul_gf264( gf.mul(y[8],v_g0), cgf.gf264_mul( cgf.index_to_gf264(offset+idx) , cgf.index_to_gf264((offset+idx)>>r1cs_dim) ) )
        values.append( cc0 )
    return gf.to_bytes(values[0])+gf.to_bytes(values[1])


def verify_proof( proof , R1CS , h_state , RS_rho = 8 , RS_shift=1<<63, verbose = 1 ) :
    if 1 == verbose : dump = print
    else : dump = _dummy

    ## process R1CS
    mat_A , mat_B , mat_C , vec_1v , witness_idx = R1CS
    m = len(mat_A)
    n = len(mat_A[0])
    pad_len = _pad_len( max(n,m) )
    dump( f"m(#rows): {m} x n: {n}, witness_idx: {witness_idx}, pad_len: {pad_len}" )

    inst_dim = _log2(witness_idx)
    r1cs_dim = _log2(pad_len)

    ## unpack proof
    rt0 = proof[0]
    rt1 = proof[1]
    _poly_len = 2*pad_len
    ldt_n_commits = fri.ldt_n_commit( _poly_len )
    ldt_commits     = proof[2:2+ldt_n_commits]
    ldt_d1poly      = proof[2+ldt_n_commits]
    ldt_open_mesgs  = proof[3+ldt_n_commits:3+ldt_n_commits+ldt_n_commits]

    open_mesgs0 = proof[3+2*ldt_n_commits]
    open_mesgs1 = proof[4+2*ldt_n_commits]

    ## recover challenges
    dump( "recover challenges" )
    h_state = H.gen( h_state , rt0 )
    chals = [ H.gen( h_state , bytes([1,i]) )[:gf.GF_BSIZE] for i in range(1,5) ]
    alpha, s1, s2, s3 = gf.from_bytes(chals[0]), gf.from_bytes(chals[1]), gf.from_bytes(chals[2]), gf.from_bytes(chals[3])
    h_state = H.gen( *chals )
    y = [ gf.from_bytes( H.gen( h_state , bytes([2,i]) )[:gf.GF_BSIZE] ) for i in range(1,10) ]
    Nq = len(open_mesgs0)
    xi, queries = fri.ldt_recover_challenges(_poly_len,h_state,ldt_commits,ldt_d1poly,Nq, RS_rho, verbose=0 )

    dump( "check if commits are opened correctly" )
    if not mt.batchverify(queries,rt0,open_mesgs0) :
        dump( "open0 fails" )
        return False
    if not mt.batchverify(queries,rt1,open_mesgs1) :
        dump( "open1 fails" )
        return False
    dump( "all passed" )

    rs_codewords = codewords_of_public_polynomials( alpha , vec_1v , mat_A , mat_B , mat_C , pad_len , RS_rho , RS_shift , verbose )

    dump( "recover first opened commit of ldt from the virtual oracle of aurora" )
    ldt_1st_mesgs = [ values_from_virtual_oracle( _idx , open_mesgs0[k][0] , open_mesgs1[k][0] , (s1,s2,s3)
                                 , y , rs_codewords , inst_dim , r1cs_dim , RS_shift ) for k,_idx in enumerate(queries) ]

    dump("verify ldt")
    ldt_r = fri.ldt_verify_proof(ldt_commits,ldt_d1poly,ldt_1st_mesgs,ldt_open_mesgs,xi,queries,RS_shift,verbose=0)
    dump( ldt_r )
    if not ldt_r : return False

    dump("all passed") 
    return True
