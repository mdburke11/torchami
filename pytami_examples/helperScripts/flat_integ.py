# this is a plain flat distribution MC integrator using a pytorch function
# with cuda

# This script will simultaneously integrate a function f_v(x) for a set of N_v
# v values (frequencies). So it will return N_v integration_results, each containing
# the individual answer, error and number of samples (N - numNans_v)

import torch
from typing import Callable


class integration_result:
    def __init__(self, sum: torch.tensor, sum2: torch.tensor, N: torch.tensor):
        self.sum = sum  #
        self.sum2 = sum2
        self.N = N

        self.ans: torch.tensor = self.sum / N
        self.avg_ans2: torch.tensor = self.sum2 / N
        self.error: torch.tensor = torch.sqrt(
            (self.avg_ans2 - self.ans**2) / self.N)


class flat_mc_integrator:
    def __init__(self, device: torch.device, N_freq, Max_batch: int = 10**6):

        self.device = device
        self.N_freq = N_freq
        self.Max_batch = Max_batch  # cutoff to do multiple batch evaluations

    def prepInput(self, length: int,
                  domain: list[list[float]]) -> torch.tensor:

        x = torch.hstack([
            torch.empty((length, 1), device=self.device).uniform_(*r)
            for r in domain
        ])
        return x

    def integrate(self, fn: Callable[[torch.tensor],
                                     torch.tensor], dim: int, N: int,
                  integration_domain: list[list[float]]) -> integration_result:

        if self.N_freq > 1:

            if dim != len(integration_domain):
                raise ValueError(
                    "Dimension and length of integration domain do not match!")

            #iters: int = N // self.Max_batch # number of extra batches to evaluate if N > batch_size

            sum: torch.tensor = torch.zeros(self.N_freq, device=self.device)
            sum2: torch.tensor = torch.zeros(self.N_freq, device=self.device)

            remaining: torch.tensor = N * torch.ones(self.N_freq,
                                                     device=self.device)

            while remaining.min() > self.Max_batch:  # take full batchs

                # while number of evals remaining exceeds max_batch size, chip
                # away with max_batch size

                # get rand vals on the domain
                x: torch.tensor = self.prepInput(self.Max_batch,
                                                 integration_domain)
                eval: torch.tensor = fn(
                    x
                )  # tensor of dim N_freq by batch_size then sum along dim=1 (rows)
                numNans = (eval != eval).sum(dim=1)
                sum += eval.nansum(dim=1)
                sum2 += (eval**2).nansum(dim=1)

                remaining = remaining - self.Max_batch + numNans

            # now only 1 batch left, so we only have the remaining
            # N - num_it * Max_batch (< Max_batch) to evaluate
            x = self.prepInput(int(remaining.min()), integration_domain)
            eval: torch.tensor = fn(x)
            numNans: torch.tensor = (eval != eval).sum(dim=1)
            sum += eval.nansum(dim=1)
            sum2 += (eval**2).nansum(dim=1)

            return integration_result(sum, sum2, N - numNans)

        else:
            # one frequency so no dim=1
            if dim != len(integration_domain):
                raise ValueError(
                    "Dimension and length of integration domain do not match!")

            #iters: int = N // self.Max_batch # number of extra batches to evaluate if N > batch_size

            sum: torch.tensor = torch.tensor(0.0, device=self.device)
            sum2: torch.tensor = torch.tensor(0.0, device=self.device)

            remaining: int = N

            while remaining > self.Max_batch:  # take full batchs

                # while number of evals remaining exceeds max_batch size, chip
                # away with max_batch size

                # get rand vals on the domain
                x: torch.tensor = self.prepInput(self.Max_batch,
                                                 integration_domain)
                eval: torch.tensor = fn(x)
                numNans = int(
                    (eval !=
                     eval).sum())  # TODO: should this be left a tensor?
                sum += eval.nansum(
                )  # TODO use torch.nansum - drops nans from batch - But then N changes!
                sum2 += (eval**2).nansum()

                remaining = remaining - self.Max_batch + numNans

            # now only 1 batch left, so we only have the remaining
            # N - num_it * Max_batch (< Max_batch) to evaluate
            x = self.prepInput(remaining, integration_domain)
            eval: torch.tensor = fn(x)
            numNans: int = int((eval != eval).sum())
            sum += eval.nansum()
            sum2 += (eval**2).nansum(
            )  # TODO use torch.nansum - drops nans from batch - But then N changes!

            return integration_result(sum, sum2, N - numNans)
