import time

import numpy as np
import scipy.sparse as sparse
from scipy.sparse import csr_matrix

from states_filtration.gauss import gauss_a_at_site
from states_filtration.gauss import get_g_squared
from states_filtration.gauss import get_g_plus

import constants as const
from tools_custom import tools_build_opers as tools
from scipy.sparse import vstack
from tools_custom.gram_schmidt import ortho_norm_vecs as gs

from scipy.linalg import null_space



N = const.STATES_NUMBER
n_sites = const.PLAQUETTE_SITES_NUMBER


def gz_filtered_states() -> list:
    zero_indices = [set(np.where(gauss_a_at_site(a=3, s_num=s_num).diagonal() == 0)[0]) for s_num in range(1, n_sites+1)]
    common_elements = zero_indices[0]
    for _set in zero_indices[1:]:
        common_elements = common_elements.intersection(_set)

    common_elements = sorted(list(common_elements))

    dimension = N ** n_sites
    gz_filtered_states: list[sparse.csr_matrix] = tools.initialize_sparse_vectors(indices=common_elements, dimension=dimension)

    return gz_filtered_states


def projector() -> csr_matrix:
    subspace_basis = vstack(gz_filtered_states())
    return subspace_basis


def g2_on_subspace() -> csr_matrix:
    g2 = get_g_squared()
    proj = projector()

    g2_subspace = proj.dot(g2.dot(proj.H))

    return g2_subspace


def g_plus_on_subspace() -> csr_matrix:
    g_plus = get_g_plus()
    proj = projector()

    g_plus_subspace = proj.dot(g_plus.dot(proj.H))
    return g_plus_subspace


# def get_phys_states() -> list: # GOT ONLY TWO ZERO EIGENVALUES -> SO IT WORKS!
#     eigenvalues, eigenvectors = np.linalg.eig(g2_on_subspace().toarray())
#
#     eigen_pairs = [(eigenvalues[i], eigenvectors[i, :]) for i in range(len(eigenvalues))]
#
#     tolerance = 1e-10
#     phys_states = []
#     for pair in eigen_pairs:
#         if pair[0] < tolerance:
#             phys_states.append(pair[1])
#
#     return gs(phys_states)


# def get_phys_states() -> list: # GOT ONLY TWO ZERO EIGENVALUES -> SO IT WORKS!
#     eigenvalues, eigenvectors = np.linalg.eig(g_plus_on_subspace().toarray())
#
#     eigen_pairs = [(eigenvalues[i], eigenvectors[i, :]) for i in range(len(eigenvalues))]
#
#     tolerance = 1e-10
#     phys_states = []
#     for pair in eigen_pairs:
#
#         if pair[0] < tolerance:
#             phys_states.append(pair[1])
#
#     return gs(phys_states)


def get_phys_states() -> list:
    g_plus = get_g_plus()
    P = projector()

    # find the null space of the matrix g_plus
    combination_matrix = g_plus @ P.T

    null_space_g_plus = null_space(combination_matrix.toarray())
    return list(null_space_g_plus.T)


def projector_phys() -> np.ndarray:
    return np.vstack(get_phys_states())


