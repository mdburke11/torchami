import torch
import pytami

class ext_vars:
    def __init__(self, beta: float, mu: complex, k: list[float]):
        self.beta = beta
        self.mu = mu
        self.k = k

    def get_ext_vars(self):
        return [self.beta, self.mu, self.k]


def epsilon_2D(
    k: torch.tensor
) -> torch.tensor:  # must return column vector of energies (use .unsqueeze(1))
    #return -2 * (torch.cos(k[:,0]) + torch.cos(k[:,1])).unsqueeze(1)
    return -2 * torch.cos(k).sum(dim=1, keepdim=True)
