# this is a plain flat distribution MC integrator using a pytorch function 
# with cuda 

import torch
from typing import Callable
import time

class integration_result:

    def __init__(self, sum: torch.tensor, sum2: torch.tensor, N: int):
        self.sum = sum
        self.sum2 = sum2
        self.N = N

        self.ans: torch.tensor = sum / N # TODO need to divide out by 1/(2pi)**dim ?
        self.avg_ans2: torch.tensor = sum2 / N
        self.error: torch.tensor = torch.sqrt((self.avg_ans2 - self.ans**2) / N)

class flat_mc_integrator:

    def __init__(self, device: torch.device, dtype=torch.float32, Max_batch: int = 10**6):

        self.device = device
        self.Max_batch = Max_batch # cutoff to do multiple batch evaluations

    def prepInput(self, length:int, domain: list[list[float]]) -> torch.tensor:
    
        x = torch.hstack([torch.empty([length, 1], device=self.device).uniform_(*r) for r in domain])
        return x

    def integrate(self, fn: Callable[[torch.tensor], torch.tensor], dim: int, N: int, 
                 integration_domain: list[list[float]]) -> integration_result:

        if dim != len(integration_domain):
            raise ValueError("Dimension and length of integration domain do not match!")
        
        iters: int = N // self.Max_batch # number of extra batches to evaluate if N > batch_size

        sum: torch.tensor = torch.tensor(0.0, device=self.device)
        sum2: torch.tensor = torch.tensor(0.0, device=self.device)

        for i in range(iters): #while num_it - iter > 1:
            # while number of evals remaining exceeds max_batch size, chip 
            # away with max_batch size

            # get rand vals on the domain
            t1 = time.time()
            x: torch.tensor = self.prepInput(self.Max_batch, integration_domain)
            t2 = time.time()
            eval: torch.tensor = fn(x)
            t3 = time.time()
            sum += eval.sum()
            sum2 += (eval**2).sum()
            print(f"PrepInput: {(t2 - t1) * 10**9} ns")
            print(f"fn: {(t3 - t2) * 10**9} ns")



        # now only 1 batch left, so we only have the remaining 
        # N - num_it * Max_batch (< Max_batch) to evaluate
        remaining: int = N - iters * self.Max_batch
        x = self.prepInput(remaining, integration_domain)

        eval: torch.tensor = fn(x)
        sum += eval.sum()
        sum2 += (eval**2).sum()

        return integration_result(sum, sum2, N)
