

from Crypto.Hash import SHA3_512 as HH

from Crypto.Random import get_random_bytes

# global prng state
class _prng(object):
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(_prng,cls).__new__(cls)
        return cls._instance

    def _update(self):
        self.buff = HH.new(self.buff).digest()

    def set_seed(self,seed=bytes(64)):
        self.used = 32
        self.buff = seed

    def gen(self, nbytes):
        r = bytearray()
        if self.used < 32 :
            ready = 32 - _prng_used
            if ready >= nbytes : ready = nbytes
            r.extend( self.buff[self.used:self.used+ready] )
            self.used += ready
            nbytes     -= ready
        while nbytes >= 32 :
            self._update()
            r.extend( self.buff )
            nbytes -= 32
        if nbytes :
            self._update()
            r.extend( self.buff[:nbytes] )
            _prng_used = nbytes
        return r


_debug = 0

#def randombytes( bsize ) :
#    return get_random_bytes( bsize )

if _debug :
    rng = _prng()
    rng.set_seed(bytes([0]*64))
    randombytes = lambda n : rng.gen(n)
else      :
    randombytes = get_random_bytes
