#include "integration.hpp"

// right now this only does real valued functions, so all things are evaluated
// twice but it should be an easy fix to make it for complex.

at::Tensor flat_mc_integrator::prepInput(int length, integ_domain &domain) {

  std::vector<at::Tensor> input_cols;
  at::TensorOptions opts =
      at::TensorOptions().dtype(at::kDouble).device(options.device());

  for (auto r : domain) {
    at::Tensor col = at::empty({length, 1}, opts);
    input_cols.push_back(at::uniform(col, r[0], r[1]));
  }
  return at::hstack(input_cols);
}

integration_result
flat_mc_integrator::integrate(std::function<at::Tensor(at::Tensor&)> fn, int dim,
                              int N, integ_domain &domain) {

  at::TensorOptions opts =
      at::TensorOptions().dtype(at::kComplexDouble).device(options.device());
  at::TensorOptions intopts =
      at::TensorOptions().dtype(at::kInt).device(options.device());

  try {
    if (domain.size() != dim) {
      throw std::invalid_argument(
          "Dimension and length of integration domain do not match!");
    }
  } catch (std::invalid_argument &e) {
    std::cerr << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }

  at::Tensor sum = at::zeros({1, N_freq}, opts);
  at::Tensor sum_2 = at::zeros({1, N_freq}, opts);

  at::Tensor remaining = N * at::ones({1, N_freq}, intopts);

  while (remaining.min().item<int>() > max_batch) {

    // while number of evals remaining exceeds max_batch size, chip
    // away with max_batch size

    // get rand vals on the domain
    at::Tensor x = this->prepInput(max_batch, domain);
    at::Tensor eval = fn(x);
    at::Tensor numNans = (eval != eval).sum(/*dim*/ 1);
    sum += eval.nansum(/*dim*/ 1);
    sum_2 += (at::pow(eval, 2)).nansum(/*dim*/ 1);

    remaining = remaining - max_batch + numNans;
  }

  // then perform one last evaluation with whats left. Note we dont re evaluate
  // for nans here, but this is accounted for.
  at::Tensor x = this->prepInput(max_batch, domain);
  at::Tensor eval = fn(x);
  at::Tensor numNans = (eval != eval).sum(/*dim*/ 1);
  sum += eval.nansum(/*dim*/ 1);
  sum_2 += (at::pow(eval, 2)).nansum(/*dim*/ 1);

  return integration_result(sum, sum_2, N - numNans);
}
