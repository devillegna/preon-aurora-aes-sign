

LAMBDA = 256

from Crypto.Hash import SHA3_256 as H

def gen( *contents ):
    return H.new(b''.join(contents)).digest()


if '__main__' == __name__ :
    print( 'hash(', b'123' , b'456' , ')->' , gen( b'123' , b'456' ) )

