# this is a plain flat distribution MC integrator using a pytorch function 
# with cuda 

import torch
from typing import Callable

class integration_result:

    def __init__(self, sum: torch.tensor, sum2: torch.tensor, N: int):
        self.sum = sum
        self.sum2 = sum2
        self.N = N

        self.ans: torch.tensor = sum / N
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
        
        #iters: int = N // self.Max_batch # number of extra batches to evaluate if N > batch_size

        sum: torch.tensor = torch.tensor(0.0, device=self.device)
        sum2: torch.tensor = torch.tensor(0.0, device=self.device)

        remaining : int = N

        while remaining > self.Max_batch: # take full batchs

            # while number of evals remaining exceeds max_batch size, chip 
            # away with max_batch size

            # get rand vals on the domain
            x: torch.tensor = self.prepInput(self.Max_batch, integration_domain)
            eval: torch.tensor = fn(x)
            numNans = int((eval != eval).sum()) # TODO: should this be left a tensor?
            sum += eval.nansum() # TODO use torch.nansum - drops nans from batch - But then N changes!
            sum2 += (eval**2).nansum()

            remaining = remaining - self.Max_batch + numNans

        # r

        # now only 1 batch left, so we only have the remaining 
        # N - num_it * Max_batch (< Max_batch) to evaluate
        x = self.prepInput(remaining, integration_domain)
        eval: torch.tensor = fn(x)
        numNans : int = int((eval != eval).sum())
        sum += eval.nansum()
        sum2 += (eval**2).nansum() # TODO use torch.nansum - drops nans from batch - But then N changes!

        return integration_result(sum, sum2, N - numNans)
