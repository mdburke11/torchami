import torch

# function returns the torch device
# if useCudaIfAvail=True, it gives the torch.cuda device if it finds one, otherwise its on the cpu
# if useCudaIfAvail=False, it defaults to the cpu
def getDevice(useCudaIfAvail=True):
    if torch.cuda.is_available() and useCudaIfAvail:
        dev = torch.device("cuda")
    else:
        dev = torch.device("cpu")
    return dev
