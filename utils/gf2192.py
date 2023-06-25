
import sys
sys.path.insert(0, '../' )

from libgf264fft import clib_wrapper as gf264

GF_EXT = 3

GF_BSIZE = 8*3

def from_bytes( bv ): return int.from_bytes(bv,byteorder='little')

def from_bytes_x2( bb ): return [ int.from_bytes(bb[:GF_BSIZE],byteorder='little'), int.from_bytes(bb[GF_BSIZE:2*GF_BSIZE],byteorder='little') ]

def to_bytes( a ): return a.to_bytes(GF_BSIZE,byteorder='little')

def from_gf264s( a0 , a1 , a2 ): return a0|(a1<<64)|(a2<<128)

def to_gf264s(a) : return a&0xffffffffffffffff,(a>>64)&0xffffffffffffffff,(a>>128)&0xffffffffffffffff

def gf2192_mul( a , b ):
    a0, a1, a2 = to_gf264s(a)
    b0, b1, b2 = to_gf264s(b)
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
