import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.fft import fft2, ifft2
# import pyfftw
# import multiprocessing
# from pyfftw.interfaces.scipy_fftpack import fft2, ifft2

# Author: Elvis do A. Soares
# Github: @elvissoares
# Date: 2020-08-16
# Updated: 2022-03-25

# pyfftw.config.NUM_THREADS = multiprocessing.cpu_count()
# print('Number of cpu cores:',multiprocessing.cpu_count())

"""
 The python script to solve the Cahn-Hilliard equation using
 an implicit pseudospectral algorithm
"""

def run(c, nsteps=1000, dt=0.1, dx=1.0, M=1.0, kappa=0.5, rho=2.0, alpha=0.0, beta=1.0):

    N = c.shape[0]
    c_hat = np.empty((N,N), dtype=np.complex64)
    dfdc_hat = np.empty((N,N), dtype=np.complex64)

    L = N * dx

    print('c0 = ',c.sum()*dx**2/L**2)

    kx = ky = np.fft.fftfreq(N, d=dx)*2*np.pi
    K = np.array(np.meshgrid(kx , ky ,indexing ='ij'), dtype=np.float32)
    K2 = np.sum(K*K,axis=0, dtype=np.float32)

    # The anti-aliasing factor
    kmax_dealias = kx.max()*2.0/3.0 # The Nyquist mode
    dealias = np.array((np.abs(K[0]) < kmax_dealias )*(np.abs(K[1]) < kmax_dealias ),dtype =bool)

    """
     The interfacial free energy density f(c) = Wc^2(1-c)^2
    """
    def finterf(c_hat):
        return kappa*ifft2(K2*c_hat**2).real

    """
     The bulk free energy density f(c) = W*c^2(1-c)^2
    """
    def fbulk(c):
        return rho * (c - alpha) ** 2 * (c - beta) ** 2

    """
     The derivative of bulk free energy density f(c) = Wc^2(1-c)^2
    """
    def dfdc(c):
        return 2 * rho * (c - alpha) * (c - beta) * (2 * c - (alpha + beta))

    def free_energy(c, c_hat):
        c_x = ifft2(c_hat * 1j * K[0]).real
        c_x_hat = fft2(c_x)
        from scipy.signal import fftconvolve
        conv = fftconvolve(c_x_hat, c_x_hat, mode='same')
        conv_ = ifft2(conv).real * dx
        other = dx * c_x**2
        print(c_x.shape)
        print(c_x_hat.shape)
        print(conv.shape)
        print(conv_.shape)
        print(conv_[:5, :5])
        print(other[:5, :5])
        raw_input('stopped')
        c_y = ifft2(c_hat * 1j * K[1]).real
        return (kappa * (c_x**2 + c_y**2) / 2. + fbulk(c)).sum() * dx**2

    c_hat[:] = fft2(c)

    c_old = c.copy()

    free_energies = []
    free_energies.append(free_energy(c, c_hat))
    for i in range(nsteps):
        dfdc_hat[:] = fft2(dfdc(c_old)) # the FT of the derivative
        dfdc_hat *= dealias # dealising
        c_hat[:] = (c_hat-dt*K2*M*dfdc_hat)/(1+dt*M*kappa*K2**2) # updating in time
        c_old[:] = c
        c = ifft2(c_hat).real # inverse fourier transform
        free_energies.append(free_energy(c, c_hat))

    print('relative_error = ',np.abs(c_old.sum()-c.sum())/c.sum())

    return c, free_energies

def plot(c, alpha, beta):
    plt.imshow(c,cmap='RdBu_r', vmin=alpha, vmax=beta)
    plt.savefig('cahn-hilliard.1f.png')
    plt.show()

def plot_free_energy(fs, nsteps, dt):
    import pandas
    dd = pandas.read_csv('moose_psu_1a_IA.csv')
    print(dd.columns)
    plt.loglog(np.arange(nsteps + 1) * dt, fs)
    plt.loglog(dd.time, dd.f_density)
    plt.show()

def initial_conc(x, y):
    return 0.5 + 0.01 * (np.cos(0.105 * x) * np.cos(0.11 * y) + (np.cos(0.13 * x) * np.cos(0.087 * y))**2 \
                         + np.cos(0.025 * x - 0.15 * y) * np.cos(0.07 * x - 0.02 * y))

if __name__ == '__main__':
    dt = 0.01
    N = 512
    Lx = 200.0
    dx  = Lx / N
    # xx = np.linspace(dx / 2.0, Lx - dx / 2.0, N)
    xx = np.linspace(0, Lx - dx, N)
    x, y = np.meshgrid(xx, xx)
    alpha = 0.3
    beta = 0.7
    c = initial_conc(x, y)
    nsteps = 10000
    # rng = np.random.default_rng(12345) # the seed of random numbers generator
    # noise = 0.1
    # c = c0 + noise * rng.standard_normal([256, 256])

    c, fs = run(c, nsteps=nsteps, dx=dx, dt=dt, alpha=0.3, beta=0.7, rho=5.0, kappa=2.0, M=5)

    # plot(c, alpha, beta)
    plot_free_energy(fs, nsteps, dt)
