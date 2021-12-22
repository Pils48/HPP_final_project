import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as si
from mpi4py import MPI


H_BAR = 1.05e-34
ELECTRON_CHARGE = 1.6e-19
K_B = 1.38e-23


class FitBessel:
    """
    This class unite important parameters for fitting by Bessel function
    Fields:
        Rn -
        Tc -
        delta_0 -
    """
    def __init__(self, Rn, Tc, delta_0):
        self.Rn = Rn
        self.Tc = Tc
        self.delta_0 = delta_0

    def Ic_ab(self, T):  # с.72 Бароне
        delta = 1.74 * self.delta_0 * np.sqrt(1. - T / self.Tc)  # [eV]

        Ch = np.pi * delta / (2 * self.Rn)
        tg = np.tanh(ELECTRON_CHARGE * delta / (2 * K_B * T))
        return Ch * tg

    def T_b(self, phi, I, T):
        Ic = self.Ic_ab(T)
        alpha = I / Ic
        gamma = (H_BAR * Ic) / (ELECTRON_CHARGE * K_B * T)

        E = np.exp(-gamma * alpha * phi / 2)
        bessel = np.i0(gamma * np.sin(phi / 2))

        return E * bessel

    def V_full(self, I, T):
        Ic = self.Ic_ab(T)
        alpha = I / Ic
        gamma = (H_BAR * Ic) / (ELECTRON_CHARGE * K_B * T)

        Integral_full = si.quad(self.T_b, 0, 2 * np.pi, args=(I, T))[0]
        P1 = 2. * self.Rn * Ic / gamma
        P2 = 1. - np.exp(-np.pi * gamma * alpha)
        P3 = 1. / Integral_full
        return P1 * P2 * P3 / 1000


def read_IV(address_folder, filename, V_inv=1, last_line=False, example=False):
    """
    Reads data from .dat file

    Args:
        address_folder : full address to a folder where the file is
        filename : full file name
        V_inv = 1 : a parameter, that flips the stress axis. if = 1 -> doesn't flip, if = -1 -> flips
        last_line = False : if True prints the last line of the file
        example = False : if True prints the 1st line of file

    Returns:
        Vocabulary: Vocabulary with sections 'Current' and 'Voltage'

    """

    address = f"{address_folder}/{filename}"
    print('Read file: ', filename)

    with open(address, 'rb') as file:
        lines = file.readlines()

    info = lines[-1:]
    if last_line:
        print('Last line: ', info[0])

    # lines = lines[:10000]
    data = lines[:-1]
    L = len(data)

    Current = []
    Voltage = []

    N = 0
    offset = 0

    for i in range(0, L):
        Measurment = (str(data[i]).replace("E", "e")).split('\\')
        Current.append(float((Measurment[0].replace("b", ""))[1:]) * 1e6)
        Voltage.append(float(Measurment[1].replace("t", "")) * 1e6 * V_inv)
        if abs(Current[i]) < 4:
            N = N + 1
            offset += Voltage[i]

    offset = offset / N
    if example:
        print(Current[1], 'mkA_______', Voltage[1], 'mkV')

    Data = {'Current': Current, 'Voltage': Voltage}

    return Data


def Split(Data):
    """
    Splits data for positive and negative sign of current

    Args:
        Data: Vocabulary with sections 'Current' and 'Voltage'

    Returns:
        I_p, V_p, I_m, V_m: 4 arrays for separete current and voltage positive and negative respectively

    """

    Current = Data['Current']
    Voltage = Data['Voltage']

    L = len(Current)
    I_p = []
    I_m = []
    V_p = []
    V_m = []

    for i in range(0, L):
        if Current[i] >= 0:
            I_p.append(Current[i])
            V_p.append(Voltage[i])
        else:
            I_m.append(Current[i] * (-1))
            V_m.append(Voltage[i] * (-1))

    return I_p, V_p, I_m, V_m


def drow_IV(Data, data_slice=0, split=True):
    """
    Drows IV plot

    Args:
        Data: Vocabulary with sections 'Current' and 'Voltage'

              or if split = True: an array with I_p, V_p, I_m, V_m.
        data_slise: X-axis starts

    Returns:
        None

    """

    if not (split):
        I = Data['Current']
        V = Data['Voltage']
        ax.plot(I, V, '.', color="red")

    I_p = []
    V_p = []
    I_m = []
    V_m = []

    for i in range(0, len(Data[0])):
        if Data[0][i] > data_slice:
            I_p.append(Data[0][i])
            V_p.append(Data[1][i])

    for i in range(0, len(Data[2])):
        if Data[2][i] > data_slice:
            I_m.append(Data[2][i])
            V_m.append(Data[3][i])

    ax.plot(I_p, V_p, linewidth=0.3, color="red", label='Positive branch')
    ax.plot(I_m, V_m, linewidth=0.3, color="green", label='Negative branch')
    ax.set_xlabel('I, $\mu$A')
    ax.set_ylabel('V, $\mu$V')
    ax.legend(shadow=True, fontsize=15)


def enumeration(s, f, d_space, T_space, R_space):
    error_min = np.inf
    delta_0_min, T_guess_min, Rn_guess_min = 0, 0, 0
    current = np.linspace(1, 9, num=100)

    for d in range(s, f):
        for T in range(s, f):
            for R in range(s, f):
                error = 0
                fit_obj = FitBessel(Rn=R_space[R], Tc=8.5, delta_0=d_space[d])
                Voltage = []
    #
                for i in range(len(current)):
                    YTR = fit_obj.V_full(current[i] * 1e-6, T_space[T]) - 0.001 * Cur[i] * 1e-6 - 0.05 * 1e-6
                    v = YTR * 1e6
                    Voltage.append(v)
                    if not np.isnan(v):
                        error += abs(correct_data[0, int(current[i] // 0.5)] - v)
                    else:
                        error += 2
                        # break
                if error < error_min and error != 0:
                    error_min = error
                    delta_0_min = d_space[d]
                    T_guess_min = T_space[T]
                    Rn_guess_min = R_space[R]

    return np.array([error_min, delta_0_min, T_guess_min, Rn_guess_min])


def search_parameters(corrected_data, **kwargs):
    resolution = kwargs["delta_resolution"]

    start_t = MPI.Wtime()
    resolution_per_rank = resolution // size

    starting_index_rank = 0
    finishing_index_rank = 0

    if rank != size:
        starting_index_rank = rank * resolution_per_rank
        finishing_index_rank = (rank + 1) * resolution_per_rank
    else:
        starting_index_rank = rank * resolution_per_rank
        finishing_index_rank = resolution_per_rank

    delta_0_space = np.linspace(kwargs["delta_s"], kwargs["delta_f"], kwargs["delta_resolution"])
    T_guess_space = np.linspace(kwargs["T_s"], kwargs["T_f"], kwargs["T_resolution"])
    Rn_guess_space = np.linspace(kwargs["Rn_s"], kwargs["Rn_f"], kwargs["Rn_resolution"])

    results = enumeration(starting_index_rank, finishing_index_rank, delta_0_space, T_guess_space, Rn_guess_space)

    recvbuf = None
    if rank == 0:
        recvbuf = np.empty([size, 4], dtype=np.float64)
    comm.Gather(results, recvbuf, root=0)

    if rank == 0:
        min_error = min(recvbuf, key=lambda arr: arr[0])
        print(min_error)
        print(f"Time elapsed {MPI.Wtime() - start_t} s")
        return min_error[1:]
    return [0, 0, 0]


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()


# Files = ['2D IV 171521.dat', '2D IV 482420.dat']
Files = ['test.dat']
address_folder = r"/home/nikita/dev/HPP/final project"
name_file = Files[0]
Data = read_IV(address_folder, name_file, -1)

delta_0 = 2.53e-3  # [eV]
T_guess = 5.63
Rn_guess = 546
Cur = np.linspace(1, 9, num=100)
correct_data = np.array([Data["Current"][:], Data["Voltage"][:]])
delta_0_1, T_guess_1, Rn_guess_1 = search_parameters(correct_data, delta_s=2.53e-3, delta_f=2.55e-3, delta_resolution=10,
                                               T_s=5.4, T_f=5.7, T_resolution=10,
                                               Rn_s=500, Rn_f=600, Rn_resolution=10)

if rank == 0:
    fig, ax = plt.subplots(figsize=(7, 10))
    ax.set_xlabel('U, $\mu V$')
    ax.set_ylabel('I, $\mu A$')
    plt.grid()

    F1 = FitBessel(Rn=Rn_guess, Tc=8.5, delta_0=delta_0)
    Voltage = []
    for i in range(len(Cur)):
        YTR = F1.V_full(Cur[i] * 1e-6, T_guess) - 0.001 * Cur[i] * 1e-6 - 0.05 * 1e-6  # - 1*Cur[i]*1e-6*np.sqrt(1-(9/Cur[i])**2)
        Voltage.append(YTR * 1e6)
    ax.plot(Cur, Voltage, label="Fit graph")

    Ic_guess = 5

    ax = drow_IV(Split(Data))
    plt.show()

