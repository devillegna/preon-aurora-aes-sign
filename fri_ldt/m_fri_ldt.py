
import sys
sys.path.insert(0, '../' )

# defining gf
from utils import gf2192 as gf

# the commit scheme
from commit import merkeltreecommit as mt

# hash for challenges
from utils import hash as H


def _dummy( *params ) : return None



####################################################

def ldt_n_commit( f_length ):
    i = 0
    while( f_length > 2 ):
        i = i+1
        f_length = f_length//2
    return i-1

###################################


def ldt_commit_phase( vi , poly_len , h_state , RS_rho=8 , verbose = 1 ):
    if 1 == verbose : dump = print
    else : dump = _dummy

    assert 0==(poly_len&(poly_len-1)) , 'poly_len is not a power of 2'
    assert len(vi)==poly_len*RS_rho   , 'RS_rho * poly_len != len(vi)'
    offset = 1<<63
    dump( f"rho: {RS_rho}, offset: {hex(offset)}" )

    commits = []
    mktrees = []

    dump( f"|original poly| : {poly_len}")
    dump( "\n### commit phase ###" )

    dump( f"Assume the evaluate on [{hex(offset)},|{poly_len}|x{RS_rho}) has been committed." )
    dump( "|v0|:" , len(vi) )
    i = 0
    while( 2 < poly_len ):
        dump( f"iteration {i}: update : [{poly_len}] -> [{poly_len//2}]:" )
        xi = gf.from_bytes( H.gen( h_state , bytes([3+i,1]) )[:gf.GF_BSIZE] )
        dump( f"derive new challenge xi <- H(h_state||{3+i}||1) : {hex(xi)}" )
        dump( f"deriving new polynomial of length [{poly_len//2}]" )

        #dump( f"v{i-1}: " , list(map(hex,vi)) )
        dump( f"ibtfy_1( |v{i}|:{len(vi)} , {hex(offset)}) ->" )
        vi = gf.ibtfy_1( vi , offset )
        #dump( f"v{i}: " , list(map(hex,vi)) )
        vi_e = vi[::2]
        vi_o = vi[1::2]
        vi = [ vi_e[j]^gf.mul(vi_o[j],xi) for j in range(len(vi_e)) ]
        dump( f"xi:{hex(xi)} * |vi|:{len(vi)} ->" )
        #dump( f"v{i}: " , list(map(hex,vi)) )
        dump( f"evaluated values generated by 1 stage of ibtfy  [{len(vi)}]" )
        offset = offset >> 1
        poly_len = poly_len//2
        i = i+1
        if poly_len <= 2 : break

        mesg = [ gf.to_bytes(vi[j]) + gf.to_bytes(vi[j+1]) for j in range(0,len(vi),2) ]
        root , randomness , tree = mt.commit( mesg )
        mktrees.append( (root,mesg,randomness,tree) )
        commits.append( root )
        dump( f"commit evaluated valuse. --> commits[{len(commits)-1}] <- |mesg|: {len(mesg)}" )
        dump( f"|commits| = {len(commits)}" )

        h_state = H.gen( h_state , gf.to_bytes(xi) , root )
        dump( f"update h_state <- H( h_state|| xi || commit ): {h_state}" )
    #cc = gf.ibtfy_1( vi , offset )     # use btfy_1 for debug only.
    #dump( "ibtfy_1:" , hex(offset) , [hex(e) for e in cc] )
    cc = gf.ifft( vi[:2] , 1 , offset )   # will get the same poly no matter applying ibtfy_1 to whatever pairs. 
    dump( "cc:" , [hex(i) for i in cc ] )
    dump( f"open deg 1 poly: {hex(cc[0])} + x* {hex(cc[1])}" )
    d1poly = gf.to_bytes(cc[0]) + gf.to_bytes(cc[1])
    h_state = H.gen( gf.to_bytes(xi) , d1poly )
    dump( f"update h_state <- H( xi || c0 || c1 ): {h_state}" )
    return commits , d1poly , mktrees , h_state



def ldt_query_phase( f_length , mktrees, h_state , Nq , RS_rho=8 , verbose = 1 ):
    if 1 == verbose : dump = print
    else : dump = _dummy
    assert Nq < 256 , "need to modify hash inputs if supporting >= 256 queries."

    dump( "\n### query phase ###" )
    dump( f"queries = [ H.gen(h_state , {3+ldt_n_commit(f_length)+1} , j )  for j in range(1,{Nq+1})]")
    queries = [ H.gen(h_state,bytes( [ 3+ldt_n_commit(f_length)+1 , j ] ))[:4] for j in range(1,Nq+1) ]   # use 32 bits of hash results only
    idx_mask = (RS_rho*f_length//2)-1
    queries = [ int.from_bytes(e,'little')&idx_mask for e in queries ]
    _queries = list(queries)
    dump( f"Queries: [{len(queries)}], {queries}" )

    # no need to open valuse of f_0
    queries = [ q//2 for q in queries ]

    open_mesgs = []
    j = 0
    for root , all_mesg, randomness, tree in mktrees :
        dump( f"open iteration: {j}" )
        open_mesgs.append( mt.batchopen(queries,all_mesg,randomness,tree) )
        dump( f"proof len:[{len(open_mesgs)}] : auth path len:{len(open_mesgs[-1][0])}" )
        queries = [ q//2 for q in queries ]
        j = j+1
    return open_mesgs , _queries


def ldt_gen_proof( f0 , h_state , Nq = 26 , RS_rho = 8 , verbose = 1 ):
    if 1 == verbose : dump = print
    else : dump = _dummy

    v0 = gf.fft( f0 , RS_rho , 1<<63 )
    dump( f"do a redundent commit for v_f0 here for checking correctness" )
    mesg0 = [ gf.to_bytes(v0[j]) + gf.to_bytes(v0[j+1]) for j in range(0,len(v0),2) ]
    rt0 , rd0 , tree0 = mt.commit( mesg0 )   # first commit

    commits , d1poly , mktrees , h_state = ldt_commit_phase( v0 , len(f0) , h_state , RS_rho , verbose )
    open_mesgs , queries = ldt_query_phase( len(f0) , mktrees , h_state , Nq , RS_rho , verbose )

    proof = [rt0]
    proof.extend( commits )
    proof.append( d1poly )
    proof.extend( open_mesgs )
    proof.append( mt.batchopen(queries,mesg0,rd0,tree0) )  # opened messages of first commit
    return proof


##########################################


def ldt_recover_challenges( _poly_len , h_state , commits , d1poly , Nq , RS_rho = 8 , verbose = 1 ):
    if 1 == verbose : dump = print
    else : dump = _dummy

    dump( "######## recovery hash state and challenges ########" )
    poly_len = _poly_len
    i = 0
    xi = []
    while( 2 < poly_len ):
        dump( f"iteration {i}: [{poly_len}] -> [{poly_len//2}]:" )
        xi.append( gf.from_bytes( H.gen( h_state , bytes([3+i,1]) )[:gf.GF_BSIZE] ) )
        dump( f"derive new challenge xi <- H(h_state||{3+i}||1) : {hex(xi[-1])}" )
        dump( f"new polynomial length [{poly_len//2}]" )
        poly_len = poly_len//2
        if poly_len <= 2 : break
        dump( f"mt.root = commits[{i}] = " ,  commits[i] )
        h_state = H.gen( h_state , gf.to_bytes(xi[i]) , commits[i] )
        dump( f"update h_state <- H( h_state|| xi || commit ): {h_state}" )
        i = i+1
    h_state = H.gen( gf.to_bytes(xi[-1]) , d1poly )
    dump( f"update h_state <- H( xi || deg1poly ): {h_state}" )

    dump( "\n### query phase ###" )
    dump( f"queries = [ H.gen(h_state , {3+i+1}=={3+ldt_n_commit(_poly_len)+1} , j )  for j in range(1,{Nq+1})]")
    queries = [ H.gen(h_state,bytes( [ 3+i+1 , j ] ))[:4] for j in range(1,Nq+1) ]
    idx_mask = (RS_rho*_poly_len//2)-1
    queries = [ int.from_bytes(e,'little')&idx_mask for e in queries ]
    dump( f"Queries: [{len(queries)}], {queries}" )
    return xi , queries



def ldt_verify_proof( commits , d1poly , first_mesgs , open_mesgs , xi , queries , Nq , verbose = 1 ):
    if 1 == verbose : dump = print
    else : dump = _dummy

    dump( "#### check linear relations and opened commit ######" )
    offset = 1<<63
    j = 0
    # check first_mesgs
    if True :
        # check linear relations
        dump( f"check linear relations:" )
        mesg      = first_mesgs   # [ path[0] for path in first_mesgs ]
        next_mesg = [ path[0] for path in open_mesgs[0] ]
        verify_j  = [ _check_linear_relation(mesg[i],next_mesg[i],queries[i],xi[j],offset) for i in range(Nq) ]
        dump( f"check linear relations:" , all(verify_j) )
        if not all(verify_j) : return False
        queries = [ q//2 for q in queries ]
        offset >>= 1
        j = j+1

    for idx,auths in enumerate(open_mesgs) :
        dump( f"open iteration: {j}" )
        dump( f"auths[{len(auths[0])}]: Nbyte: ", sum( map( len,auths[0]) ) )
        if not mt.batchverify( queries , commits[j-1] , auths ) :
            dump("batchverify() fails")
            return False
        else : dump("oepned mesgs are verified.")

        # check linear relations
        mesg = [ path[0] for path in auths ]
        if idx == len(open_mesgs)-1 : break
        dump( f"check linear relations [{idx}]:" )
        next_mesg = [ path[0] for path in open_mesgs[idx+1] ]
        verify_j  = [ _check_linear_relation(mesg[i],next_mesg[i],queries[i],xi[j],offset) for i in range(Nq) ]
        dump( f"check linear relations [{idx}]:" , all(verify_j) )
        if not all(verify_j) : return False
        queries = [ q//2 for q in queries ]
        offset >>= 1
        j = j+1
    # check deg 1 poly
    verify_j = [ _check_deg1poly_linear_relation(mesg[i],d1poly,queries[i],xi[-1],offset) for i in range(Nq) ]
    dump( f"check last linear relations (with the d1poly):" , all(verify_j) )
    if not all(verify_j) : return False
    return True

def _check_linear_relation( mesgj1 , mesgj0 , idx , xi , offset ) :
    new_j0 = gf.from_bytes_x2( mesgj0 )
    org_j1 = gf.from_bytes_x2( mesgj1 )
    org_j0 = gf.ibtfy_1( org_j1 , offset^(idx<<1) )
    cc1 = org_j0[0] ^ gf.mul( org_j0[1] , xi )
    return new_j0[idx&1] == cc1

def _check_deg1poly_linear_relation( mesgjm1 , d1poly , idx , xi , offset ) :
    m0 = gf.fft( gf.from_bytes_x2(d1poly) , 1 , (offset>>1)^(idx^(idx&1)) )
    return _check_linear_relation( mesgjm1 , gf.to_bytes(m0[0])+gf.to_bytes(m0[1]) , idx , xi , offset )


def ldt_verify( proof , _poly_len , h_state , Nq = 26 , RS_rho = 8 , verbose = 1 ):
    n_commits = ldt_n_commit( _poly_len )
    first_commit = proof[0]
    commits     = proof[1:1+n_commits]
    d1poly      = proof[1+n_commits]
    open_mesgs  = proof[2+n_commits:2+n_commits+n_commits]
    first_mesgs = proof[2+n_commits+n_commits]
    xi, queries = ldt_recover_challenges(_poly_len,h_state,commits,d1poly,Nq, RS_rho, verbose )

    if not mt.batchverify( queries , first_commit , first_mesgs ) : return False

    return ldt_verify_proof(commits,d1poly,[ path[0] for path in first_mesgs ],open_mesgs,xi,queries,Nq,verbose)

