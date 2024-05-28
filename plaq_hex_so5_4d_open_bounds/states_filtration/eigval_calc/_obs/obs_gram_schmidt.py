import numpy as np
import constants as const
from states_filtration import gauss
from ham_electric.ham_el import ham_el_plaq
from ham_magnetic.ham_mag import ham_mag_plaq


def compute_eigvals2(g, q=const.TRUNCATION) -> np.ndarray:
    ham_el = ham_el_plaq(g=g, q=q)
    ham_mag = ham_mag_plaq(g=g, q=q)

    phys_states = gauss.get_phys_states(q=q)

    v1, v2 = phys_states[0], phys_states[1]

    w1 = v1 / np.linalg.norm(v1)

    w2 = v2 - np.dot(np.conj(v1), v2) * v1 / (np.linalg.norm(v1)**2)
    w2 = w2 / np.linalg.norm(w2)

    w = [w1, w2]
    ham_full = ham_el + ham_mag
    ham_phys = [[np.dot(np.conj(wi), (ham_full @ wj)) for wj in w] for wi in w]
    ham_phys = np.array(ham_phys)

    eigenvalues = np.linalg.eig(ham_phys)[0]

    return eigenvalues


def compute_eigvals3(g, q=1) -> np.ndarray:
    ham_el = ham_el_plaq(g=g, q=q)
    ham_mag = ham_mag_plaq(g=g, q=q)

    phys_states = gauss.get_phys_states(q=q)

    v1, v2, v3 = phys_states[0], phys_states[1], phys_states[2]

    w1 = v1 / np.linalg.norm(v1)
    w2 = v2 - np.dot(np.conj(w1), v2) * w1 / (np.linalg.norm(w1) ** 2)
    w2 = w2 / np.linalg.norm(w2)

    w3 = v3 - np.dot(np.conj(w1), v3) * w1 / (np.linalg.norm(w1) ** 2) - \
    np.dot(np.conj(w2), v3) * w2 / (np.linalg.norm(w2) ** 2)
    w3 = w3 / np.linalg.norm(w3)

    w = np.array([w1, w2, w3])
    ham_mag_phys = [[np.dot(np.conj(wi), (ham_mag @ wj)) for wj in w] for wi in w]
    ham_mag_phys = np.array(ham_mag_phys)

    ham_el_phys = [[np.dot(np.conj(wi), (ham_el @ wj)) for wj in w] for wi in w]
    ham_el_phys = np.array(ham_el_phys)

    ham_phys = ham_el_phys + ham_mag_phys
    eigvals = np.sort(np.linalg.eig(ham_phys)[0])

    return eigvals