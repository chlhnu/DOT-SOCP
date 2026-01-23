import numpy as np
from .initialize import VariableState, ModelState


def norm_l2(x, h):
    return np.sqrt(h) * np.linalg.norm(x)


def initial_scaling(var: VariableState, model: ModelState, scaling_yes=True, last_level_kkt=None):
    c = model.c
    grad = model.grad

    h = 1.0 / var.phi.shape[0]
    h_mean = h ** (1.0 / (model.dim + 1))

    if last_level_kkt is None or not hasattr(var, 'E2'):
        e_scale2 = np.sqrt(2.0)
    else:

        ratio = np.sqrt(last_level_kkt[0] / (last_level_kkt[1] + 1e-10))
        lower_ratio = 0.8333

        current_e2 = getattr(var, 'E2', np.sqrt(2.0))

        if ratio < lower_ratio:
            e_scale2 = current_e2 * max(1.0 / np.sqrt(2.0), ratio / lower_ratio)
        else:
            e_scale2 = current_e2 * min(np.sqrt(2.0), max(1.0, ratio))

    if scaling_yes:

        norm_c_val = norm_l2(c, h) * np.sqrt(model.nt)
        norm_d_val = np.sqrt(2.0)

        adjust = 1.0 if model.weight is None else np.pow(10, np.mean(np.log10(model.weight + 1e-10)))
        D = np.sqrt(2.0) * np.sqrt(h_mean) * adjust
        E = D / e_scale2

        c_scale = max(1.0, norm_c_val * np.sqrt(h_mean) / adjust)
        d_scale = E * norm_d_val * np.sqrt(adjust)

        model.c = c / c_scale
        model.normc = norm_c_val / c_scale
        model.normd = norm_d_val * E / d_scale
        model.grad = D * grad

        var.phi = (1.0 / d_scale) * var.phi
        var.q = (D / d_scale) * var.q
        var.z = (E / d_scale) * var.z
        var.alpha = (1.0 / (c_scale * D)) * var.alpha
        var.beta = (1.0 / (c_scale * E)) * var.beta

    else:
        c_scale = 1.0
        d_scale = 1.0
        D = 1.0
        E = 1.0
        e_scale2 = np.sqrt(2.0)

        model.normc = norm_l2(c, h)
        model.normd = np.sqrt(2.0)


    var.cScale = c_scale
    var.dScale = d_scale
    var.D = D
    var.E = E
    var.E2 = e_scale2

    return var, model


def recover_org_var(var: VariableState):
    c_scale = var.cScale
    d_scale = var.dScale
    D = var.D
    E = var.E
    
    assert c_scale is not None and d_scale is not None and D is not None and E is not None, "Scaling factors are not set in VariableState."

    var.phi = d_scale * var.phi
    var.z = (d_scale / E) * var.z
    var.q = (d_scale / D) * var.q
    var.alpha = (c_scale * D) * var.alpha
    var.beta = (c_scale * E) * var.beta

    return var