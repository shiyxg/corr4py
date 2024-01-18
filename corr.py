import numpy as np
import os
import sys
import time
pwd = os.getcwd()
sys.path.append(pwd+'/lib')

THREADS_NUM = 8

from corrFFT import corr_fft_full
from corrFFT import corr_fft_padding_f64
from corrFFT import corr_fft_padding_f32
from corrFFT import corr_fft_valid_f64
from corrFFT import corr_fft_valid_f32

def corr_with_wavelets(data:np.array,wavelets:np.array, sp=None, ep=None,
                       THREADS_NUM=THREADS_NUM) -> np.array:
    """
    data: ne,ns,nt
    wavelets: ne,nt
    sp: valid mini shift value for corr, default: -nt+1
    ep: valid mini shift value for corr, default: nt

    return: data[i,j,:] x wavelets[i,:]*, of which shape=[ne,ns,ep-sp]
    """
    ne,ns,nt = data.shape
    ne1, nt1 = wavelets.shape
    assert ne1==ne and nt1==nt

    if sp is None:
        sp = -nt+1
    if ep is None:
        ep = nt

    assert type(sp)==type(1)
    assert type(ep)==type(1)
    
    if sp<-nt+1 or sp>nt-1 or ep>nt or ep<=-nt+1 or ep<=sp:
        raise ValueError(f'sp={sp},ep={ep}:index of [{sp},{ep}) is not valid for corr:[{-nt+1},{nt})')

    d_1d = np.reshape(data,[ne*ns*nt])
    w_1d = np.reshape(wavelets,[ne*nt])
    if data.dtype == np.float32:
        c1 = corr_fft_valid_f32(d_1d, w_1d, ne,ns,nt,sp,ep,THREADS_NUM)
    if data.dtype == np.float64:
        c1 = corr_fft_valid_f64(d_1d, w_1d, ne,ns,nt,sp,ep,THREADS_NUM)

    c_3d = np.reshape(c1, [ne,ns,ep-sp],order='C')

    return c_3d

def corr_fft(data:np.array,wavelets:np.array, sp=None, ep=None,
                       THREADS_NUM=THREADS_NUM) -> np.array:
    """
    data: ne,nt
    wavelets: ne,nt
    sp: valid mini shift value for corr, default: -nt+1
    ep: valid mini shift value for corr, default: nt

    return: data[i,:] x wavelets[i,:]*, of which shape=[ne,ep-sp]
    """

    ne, nt = data.shape
    ne1, nt1 = wavelets.shape
    assert ne1==ne and nt1==nt

    if sp is None:
        sp = -nt+1
    if ep is None:
        ep = nt

    assert type(sp)==type(1)
    assert type(ep)==type(1)
    
    if sp<-nt+1 or sp>nt-1 or ep>nt or ep<=-nt+1 or ep<=sp:
        raise ValueError(f'sp={sp},ep={ep}:index of [{sp},{ep}) is not valid for corr:[{-nt+1},{nt})')

    d_1d = np.reshape(data,[ne*nt])
    w_1d = np.reshape(wavelets,[ne*nt])
    ns=1
    if data.dtype == np.float32:
        c1 = corr_fft_valid_f32(d_1d, w_1d, ne,1,nt,sp,ep,THREADS_NUM)
    if data.dtype == np.float64:
        c1 = corr_fft_valid_f64(d_1d, w_1d, ne,1,nt,sp,ep,THREADS_NUM)

    c_3d = np.reshape(c1, [ne,ns,ep-sp],order='C')

    return c_3d
    
def benchmark(func='F64'):
    print('init matrix')
    ne,ns,nt = 1000,100,10000
    a = np.random.random([ne,ns,nt])
    b = np.random.random([ne,nt])
    b[:,0]=1
    s = time.time()
    print('corr4py start')
    if func=='F64':
        # sp = -nt+1
        # ep = nt
        # c1 = corr_fft_valid_f64(np.reshape(a,[ne*ns*nt]), np.reshape(b,[ne*nt]), ne,ns,nt,sp,ep,THREADS_NUM)
        # c1 = np.reshape(c1, [ne,ns,ep-sp])
        c1=corr_with_wavelets(a,b)
    elif func=='F32':
        a = a.astype('float32')
        b = b.astype('float32')
        # sp = -nt+1
        # ep = nt
        # c1 = corr_fft_valid_f32(np.reshape(a,[ne*ns*nt]), np.reshape(b,[ne*nt]), ne,ns,nt,sp,ep,THREADS_NUM)
        # c1 = np.reshape(c1, [ne,ns,ep-sp])
        c1=corr_with_wavelets(a,b)
    else:
        raise ValueError(f'func={func} not in [F32, F64]')
        # return 0
    e = time.time()
    print(f'Matrix:{ne} x {ns} x {nt} x {func}', )
    print('time of coor4py:',e-s,'s')

    from scipy import signal

    s = time.time()
    print('scipy start')
    if func=='F64':
        c2 = np.zeros([ne,ns,nt*2-1])
    else:
        c2 = np.zeros([ne,ns,nt*2-1],dtype='float32')
    for i in range(ne):
        for j in range(ns):
            c2[i,j,:] = signal.correlate(a[i,j,:],b[i,:])
        print(f'\r {i}|{ne}',end='')
    e = time.time()
    print('\n time of scipy',e-s,'s')

    if ne*ns*nt<10:
        print(a)
        print(c1)
        print(c2)

    print('summary difference between coor4py and scipy.signal.correlate:')
    print(np.sum(np.abs(c1-c2)))

    return 1


if __name__ =='__main__':
    print('start benchmark of F64 with {} threads'.format(THREADS_NUM))
    benchmark('F64')
    print('#'*15)
    print('start benchmark of F32 with {} threads'.format(THREADS_NUM))
    benchmark('F32')