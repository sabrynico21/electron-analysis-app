"""Latent correlation estimation (port of the `latentcor` R package)."""

from .latentcor import LatentCorResult, latentcor

__all__ = ["LatentCorResult", "latentcor"]
