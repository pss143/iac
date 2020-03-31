
library(yaml)

#' Read the network specificaton from a YAML file
#'
#' Details about how to setup the yaml file forthcoming
#' @param yaml_file Name of yaml file
#' @param verbose If true, then pools and connections are printed to the screen (default=FALSE)
#' @keywords network read
#' @export
#' @examples
#' iac.readnet('dualroutes.yaml', verbose = TRUE)
#'
iac.readnet = function(yaml_file, verbose = FALSE) {

  setslices = function(a, b, ai, bi, w, dist = 'cityblock') {
    stopifnot(any(startsWith(c('euclidean', 'cityblock'), dist)))
    new.dims = c(a, b)
    ndims = length(new.dims)
    arr = array(0, dim = new.dims)
    ai.i = ai
    bi.i = bi + length((a))
    stopifnot(identical(new.dims[ai.i], new.dims[bi.i]))
    dimtypes = 1:ndims
    dimtypes[ai.i] = 'mapped.a'
    dimtypes[bi.i] = 'mapped.b'
    dimtypes[-c(ai.i, bi.i)] = 'free'
    corres = cbind(ai.i, bi.i)
    specs = list()
    for(i in seq(length(dimtypes))) {
      specs[[i]] = seq(new.dims[i])
    }
    cxn.list = expand.grid(specs)
    distances = cbind(cxn.list[, corres[, 1]] - cxn.list[, corres[, 2]], 0)
    euclid = sqrt(rowSums(distances^2))
    cityblock = rowSums(abs(distances))
    for(i in seq(nrow(cxn.list))) {
      v = as.numeric(cxn.list[i, ])
      d = ifelse(startsWith('euclid', dist), euclid[i],
                 ifelse(startsWith('cityblock', dist), cityblock[i], NA))
      d = floor(d + 1)
      weight = ifelse(d > length(w), 0, w[d])
      arr[matrix(v, 1)] = weight
    }

    # return absoluately everything! Might be handy for debugging
    return(list(grid = cxn.list, arr = arr,
                specs= specs, corres = corres,
                euclid = euclid, city = cityblock))
  }

  connect.pools = function(network, from.name, to.name,
                           from.dims=NULL, to.dims=NULL, cxn.wts,
                           mode='set', dist = 'city', reciprocal = FALSE) {
    from = network$pools[[from.name]]
    to = network$pools[[to.name]]
    if(is.null(from.dims)) { from.dims = seq(length(from$shape)) }
    if(is.null(to.dims))   { to.dims   = seq(length(to$shape)) }

    res = setslices(from$shape, to$shape, from.dims, to.dims,
                    w = cxn.wts, dist = dist)
    fromi = rep(from$units, to$size)
    toi = rep(to$units, each = from$size)
    wtvec = as.vector(res$arr)
    for(i in 1:length(wtvec)) {
      if(mode == 'set') {
        network$wt.matrix[fromi[i], toi[i]] = wtvec[i]
        if(reciprocal) {
          network$wt.matrix[toi[i], fromi[i]] = wtvec[i]
        }
      }
      else if(mode == 'add') {
        network$wt.matrix[fromi[i], toi[i]] = wtvec[i] + network$wt.matrix[fromi[i], toi[i]]
        if(reciprocal) {
          network$wt.matrix[toi[i], fromi[i]] = wtvec[i] + network$wt.matrix[toi[i], fromi[i]]
        }
      }
    }
    return(network)
  }

  mknet = function(poolinfo) {
    pools = {}
    poolnames = names(poolinfo)
    nunits = 0
    for(i in 1:length(poolinfo)) {
      pool = list(shape = poolinfo[[i]])
      pool$size = prod(pool$shape)
      pool$units = (1:pool$size) + nunits
      pools[[poolnames[i]]] = pool
      nunits = nunits + pool$size
    }
    network = list(nunits = nunits,
                   pools = pools,
                   wt.matrix = matrix(0, nrow = nunits, ncol = nunits))
    return(network)
  }

  spec = yaml.load_file(yaml_file)
  net = mknet(spec$pools)
  for(cxn in spec$connections) {
    if(is.null(cxn$reciprocal)) { cxn$reciprocal = FALSE }
    cxn$from.dims = unlist(cxn$from.dims)
    cxn$to.dims = unlist(cxn$to.dims)
    cxn$w = unlist(cxn$w)
    if(verbose) { message(str(cxn)) }
    net = connect.pools(net, from.name = cxn$from, to.name = cxn$to,
                        from.dims = unlist(cxn$from.dims),
                        to.dims = unlist(cxn$to.dims),
                        cxn.wts = unlist(cxn$w),
                        reciprocal = cxn$reciprocal)
  }
  return(net)
}

#' Read the network parameters from a YAML file
#' Expected to specify the following (example values shown):
# alpha: 0.25     # weighting of excitatory inputs
# gamma: 0.25     # weighting of inhibitory inputs
# ext.weight: 0.4 # weighting of external input
# decay: 0.25     # speed at which activations return to resting value
# noise: 0.01     # sd of gaussian noise added to net input (0 to disable)
# min: -1.0       # minimum activation value (just leave as -1)
# max:  1.0       # maximum actiation value (just leave as +1)
# rest: 0.0       # resting activation level, no reason to change
#'
#' @param yaml_file Name of yaml file
#' @keywords network read
#' @export
#' @examples
#' iac.readparams('params.yaml')
#'
iac.readparams = function(yaml_file) {
  params = yaml.load_file(yaml_file)
  return(params)
}

#' Cycle -- run the network
#'
#' Run one cycle of the network
#' Usually you would do this after creating the network, eg net=iac.readfile('dualroutes.yaml')
#' @param state.act A vector of current unit activations (length=net$nunits)
#' @param weights Weight matrix, normally net$wt.matrix
#' @param external.act  A vector of external input (length=net$nunits)
#' @param paramlist A list(noise=x, alpha=x, gamma=x, ext.weight=x, rest=x, min=x, max=x, decay=x).Normally read this in from a file with iac.readparams()
#' @keywords network run
#' @export
#' @examples
#' state.act = iac.cycle(state.act, net$wt.matrix, external.act, params)
#' log.act = iac.log.update(log.act, state.act)
iac.cycle = function(state.act, weights, external.act, paramlist) {
  getinput = function(state.act, weights, external.act, external.weight,
                      alpha, gamma, noise) {
    nunits = length(state.act)
    weights.positive = ifelse(weights > 0, weights, 0)
    weights.negative = ifelse(weights < 0, weights, 0)
    net.external = external.act * external.weight
    net.excitation = colSums(weights.positive * state.act) * alpha
    net.inhibition = colSums(weights.negative * state.act) * gamma
    input = net.external + net.excitation + net.inhibition
    if(noise > 0) { input = input + rnorm(nunits, 0, noise) }
    return(input)
  }
  update = function(state.act, input, min, max, rest, decay) {
    input.positive = ifelse(input > 0, input, 0)
    input.negative = ifelse(input < 0, input, 0)
    delta = (max - state.act) * input.positive +
      (state.act - min) * input.negative
    delta = delta - decay * (state.act - rest)
    state.act = state.act + delta
  }
  input = getinput(state.act, weights, external.act,
                   paramlist[['ext.weight']],
                   paramlist[['alpha']], paramlist[['gamma']],
                   paramlist[['noise']])
  state.act = update(state.act, input,
                     paramlist[['min']], paramlist[['max']],
                     paramlist[['rest']], paramlist[['decay']])
  return(state.act)
}

#' Update the log of activations
#'
#' Usually you would do this after every cycle of the network
#' @param log.act A matrix cols=net$nunits, rows=activations on each cycle
#' @param state.act The activations just computed from iac.cycle
#' @keywords network run
#' @export
#' @examples
#' state.act = iac.cycle(state.act, net$wt.matrix, external.act, .4, params)
#' log.act = iac.log.update(log.act, state.act)

iac.log.update = function(current.log, new.state.act) {
  current.log = rbind(current.log, new.state.act, deparse.level = 0)
  dimnames(current.log)[[1]] = 0:(nrow(current.log) - 1)
  return(current.log)
}

