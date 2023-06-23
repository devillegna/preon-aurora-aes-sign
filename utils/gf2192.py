
import sys
sys.path.insert(0, '../' )

from libgf264fft import clib_wrapper as gf264

GF_EXT = 3

def gf2192_mul( a , b ):
    a0 = a&0xffffffffffffffff
    a1 = (a>>64)&0xffffffffffffffff
    a2 = (a>>128)&0xffffffffffffffff
    b0 = b&0xffffffffffffffff
    b1 = (b>>64)&0xffffffffffffffff
    b2 = (b>>128)&0xffffffffffffffff
    c0 = gf264.gf264_mul(a0,b0)
    c1 = gf264.gf264_mul(a0,b1)^gf264.gf264_mul(a1,b0)
    c2 = gf264.gf264_mul(a0,b2)^gf264.gf264_mul(a1,b1)^gf264.gf264_mul(a2,b0)
    c3 = gf264.gf264_mul(a1,b2)^gf264.gf264_mul(a2,b1)
    c4 = gf264.gf264_mul(a2,b2)
    c2 ^= c4
    c1 ^= c4^c3
    c0 ^= c3
    return c0|(c1<<64)|(c2<<128)

def mul( a , b ) : return gf2192_mul( a , b )
