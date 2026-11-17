"""Stable identities for independent uniform sensitivity parameters."""

from dataclasses import asdict, dataclass
import hashlib
import json


@dataclass(frozen=True)
class Parameter:
    """A native target, incident interval and shape mode (energies in eV)."""

    nuclide: str
    mt: int
    family: str
    mode: int = 0
    product: int = 0
    distribution: int = 0
    interval: tuple[float, float] = (0.0, 20e6)
    limit: float = 0.2
    shape: str = "strict"
    xs_mode: str = "linear"

    def __post_init__(self):
        if self.family not in ("xs", "energy", "angle"):
            raise ValueError("Unknown family")
        if (
            not 0 < self.limit < 1
            or not 0 <= self.interval[0] < self.interval[1]
        ):
            raise ValueError("Invalid bounds or incident interval")
        if self.shape not in ("strict", "normalized") or self.xs_mode not in (
            "linear",
            "exponential",
        ):
            raise ValueError("Unknown perturbation semantics")
        if self.mode not in ((0,) if self.family == "xs" else (1, 2)):
            raise ValueError("Unsupported mode")

    @property
    def identity(self) -> str:
        """Return stable identity including basis and bounds."""
        digest = hashlib.sha256(
            json.dumps(asdict(self), sort_keys=True).encode()
        ).hexdigest()[:12]
        return f"{self.nuclide}:MT{self.mt}:{self.family}:p{self.product}:d{self.distribution}:m{self.mode}:{digest}"

    @property
    def target(self) -> tuple:
        return (
            self.nuclide,
            self.mt,
            self.product,
            self.distribution,
            self.family,
        )

    @property
    def reaction(self) -> str:
        return f"{self.nuclide}:MT{self.mt}"

    def mask(self, grid):
        return (grid >= self.interval[0]) & (grid <= self.interval[1])

    @property
    def basis_definition(self) -> str:
        """Document the native-grid basis in exported registries."""
        if self.family == "xs":
            return (
                "Native incident-interval indicator; sigma*(1+delta)"
                if self.xs_mode == "linear"
                else "Native incident-interval indicator; sigma*exp(delta)"
            )
        function = (
            ("z=Eout/max(Eout)" if self.mode == 1 else "z², z=Eout/max(Eout)")
            if self.family == "energy"
            else ("P1(mu)=mu" if self.mode == 1 else "P2(mu)=(3mu²-1)/2")
        )
        return (
            function
            + "; "
            + self.shape
            + " normalized shape; two modes share budget"
        )
