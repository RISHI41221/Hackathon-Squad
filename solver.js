"use strict";

const fs = require("fs");

function makeScanner(text) {
  let idx = 0;
  const n = text.length;
  return function nextInt() {
    while (idx < n) {
      const c = text.charCodeAt(idx);
      if (c > 32) {
        break;
      }
      idx += 1;
    }
    let sign = 1;
    if (text.charCodeAt(idx) === 45) {
      sign = -1;
      idx += 1;
    }
    let value = 0;
    while (idx < n) {
      const c = text.charCodeAt(idx);
      if (c < 48 || c > 57) {
        break;
      }
      value = value * 10 + (c - 48);
      idx += 1;
    }
    return value * sign;
  };
}

function parseInput(text) {
  const nextInt = makeScanner(text);
  const n = nextInt();
  const m = nextInt();
  const weights = new Float64Array(n);
  for (let i = 0; i < n; i += 1) {
    weights[i] = nextInt();
  }

  const deg = new Int32Array(n);
  const edgesU = new Int32Array(m);
  const edgesV = new Int32Array(m);

  for (let i = 0; i < m; i += 1) {
    const u = nextInt() - 1;
    const v = nextInt() - 1;
    edgesU[i] = u;
    edgesV[i] = v;
    deg[u] += 1;
    deg[v] += 1;
  }

  const start = new Int32Array(n + 1);
  for (let i = 0; i < n; i += 1) {
    start[i + 1] = start[i] + deg[i];
  }

  const adj = new Int32Array(2 * m);
  const ptr = new Int32Array(start);
  for (let i = 0; i < m; i += 1) {
    const u = edgesU[i];
    const v = edgesV[i];
    adj[ptr[u]++] = v;
    adj[ptr[v]++] = u;
  }

  return { n, m, weights, deg, start, adj, edgesU, edgesV };
}

function setLocalIndex(vertices, localIndex) {
  for (let i = 0; i < vertices.length; i += 1) {
    localIndex[vertices[i]] = i;
  }
}

function clearLocalIndex(vertices, localIndex) {
  for (let i = 0; i < vertices.length; i += 1) {
    localIndex[vertices[i]] = -1;
  }
}

function applyTreeSolution(instance, vertices, localIndex, chosenGlobal) {
  const { weights, start, adj } = instance;
  const k = vertices.length;
  setLocalIndex(vertices, localIndex);

  const parent = new Int32Array(k);
  parent.fill(-1);
  const order = new Int32Array(k);
  let head = 0;
  let tail = 0;
  const queue = new Int32Array(k);

  queue[tail++] = 0;
  parent[0] = 0;
  let orderSize = 0;

  while (head < tail) {
    const li = queue[head++];
    order[orderSize++] = li;
    const gv = vertices[li];
    for (let e = start[gv]; e < start[gv + 1]; e += 1) {
      const to = adj[e];
      const lj = localIndex[to];
      if (lj === -1 || parent[lj] !== -1) {
        continue;
      }
      parent[lj] = li;
      queue[tail++] = lj;
    }
  }

  const take = new Float64Array(k);
  const skip = new Float64Array(k);
  for (let idx = k - 1; idx >= 0; idx -= 1) {
    const li = order[idx];
    const gv = vertices[li];
    let takeValue = weights[gv];
    let skipValue = 0;
    for (let e = start[gv]; e < start[gv + 1]; e += 1) {
      const to = adj[e];
      const child = localIndex[to];
      if (child === -1 || parent[child] !== li) {
        continue;
      }
      takeValue += skip[child];
      skipValue += Math.max(take[child], skip[child]);
    }
    take[li] = takeValue;
    skip[li] = skipValue;
  }

  const chosenLocal = new Uint8Array(k);
  const stackV = new Int32Array(k);
  const stackTaken = new Uint8Array(k);
  let sp = 0;
  stackV[sp] = 0;
  stackTaken[sp] = 0;
  sp += 1;

  while (sp > 0) {
    sp -= 1;
    const li = stackV[sp];
    const parentTaken = stackTaken[sp] === 1;
    const takeSelf = !parentTaken && take[li] >= skip[li];
    if (takeSelf) {
      chosenLocal[li] = 1;
      chosenGlobal[vertices[li]] = 1;
    }
    const gv = vertices[li];
    for (let e = start[gv]; e < start[gv + 1]; e += 1) {
      const to = adj[e];
      const child = localIndex[to];
      if (child === -1 || parent[child] !== li) {
        continue;
      }
      stackV[sp] = child;
      stackTaken[sp] = takeSelf ? 1 : 0;
      sp += 1;
    }
  }

  clearLocalIndex(vertices, localIndex);
}

function bipartiteColoring(instance, vertices, localIndex) {
  const { start, adj } = instance;
  const k = vertices.length;
  setLocalIndex(vertices, localIndex);

  const color = new Int8Array(k);
  color.fill(-1);
  const queue = new Int32Array(k);

  for (let root = 0; root < k; root += 1) {
    if (color[root] !== -1) {
      continue;
    }
    let head = 0;
    let tail = 0;
    color[root] = 0;
    queue[tail++] = root;
    while (head < tail) {
      const li = queue[head++];
      const gv = vertices[li];
      for (let e = start[gv]; e < start[gv + 1]; e += 1) {
        const lj = localIndex[adj[e]];
        if (lj === -1) {
          continue;
        }
        if (color[lj] === -1) {
          color[lj] = color[li] ^ 1;
          queue[tail++] = lj;
        } else if (color[lj] === color[li]) {
          clearLocalIndex(vertices, localIndex);
          return null;
        }
      }
    }
  }

  return color;
}

class Dinic {
  constructor(n) {
    this.n = n;
    this.head = new Int32Array(n);
    this.head.fill(-1);
    this.to = [];
    this.next = [];
    this.cap = [];
    this.level = new Int32Array(n);
    this.ptr = new Int32Array(n);
  }

  addEdge(u, v, c) {
    this.to.push(v);
    this.cap.push(c);
    this.next.push(this.head[u]);
    this.head[u] = this.to.length - 1;

    this.to.push(u);
    this.cap.push(0);
    this.next.push(this.head[v]);
    this.head[v] = this.to.length - 1;
  }

  bfs(source, sink) {
    this.level.fill(-1);
    const queue = new Int32Array(this.n);
    let head = 0;
    let tail = 0;
    this.level[source] = 0;
    queue[tail++] = source;
    while (head < tail) {
      const v = queue[head++];
      for (let e = this.head[v]; e !== -1; e = this.next[e]) {
        if (this.cap[e] > 0 && this.level[this.to[e]] === -1) {
          this.level[this.to[e]] = this.level[v] + 1;
          queue[tail++] = this.to[e];
        }
      }
    }
    return this.level[sink] !== -1;
  }

  dfs(v, sink, pushed) {
    if (pushed === 0) {
      return 0;
    }
    if (v === sink) {
      return pushed;
    }
    for (let e = this.ptr[v]; e !== -1; e = this.next[e], this.ptr[v] = e) {
      if (this.cap[e] <= 0 || this.level[this.to[e]] !== this.level[v] + 1) {
        continue;
      }
      const tr = this.dfs(this.to[e], sink, Math.min(pushed, this.cap[e]));
      if (tr === 0) {
        continue;
      }
      this.cap[e] -= tr;
      this.cap[e ^ 1] += tr;
      return tr;
    }
    return 0;
  }

  maxFlow(source, sink) {
    let flow = 0;
    while (this.bfs(source, sink)) {
      for (let i = 0; i < this.n; i += 1) {
        this.ptr[i] = this.head[i];
      }
      while (true) {
        const pushed = this.dfs(source, sink, Number.MAX_SAFE_INTEGER);
        if (pushed === 0) {
          break;
        }
        flow += pushed;
      }
    }
    return flow;
  }

  reachable(source) {
    const seen = new Uint8Array(this.n);
    const queue = new Int32Array(this.n);
    let head = 0;
    let tail = 0;
    queue[tail++] = source;
    seen[source] = 1;
    while (head < tail) {
      const v = queue[head++];
      for (let e = this.head[v]; e !== -1; e = this.next[e]) {
        if (this.cap[e] > 0 && seen[this.to[e]] === 0) {
          seen[this.to[e]] = 1;
          queue[tail++] = this.to[e];
        }
      }
    }
    return seen;
  }
}

function applyBipartiteSolution(instance, vertices, localIndex, color, chosenGlobal) {
  const { weights, start, adj } = instance;
  const k = vertices.length;
  const source = k;
  const sink = k + 1;
  const dinic = new Dinic(k + 2);
  let totalWeight = 0;
  const inf = Number.MAX_SAFE_INTEGER / 8;

  for (let i = 0; i < k; i += 1) {
    const gv = vertices[i];
    totalWeight += weights[gv];
    if (color[i] === 0) {
      dinic.addEdge(source, i, weights[gv]);
    } else {
      dinic.addEdge(i, sink, weights[gv]);
    }
  }

  for (let i = 0; i < k; i += 1) {
    const gv = vertices[i];
    for (let e = start[gv]; e < start[gv + 1]; e += 1) {
      const j = localIndex[adj[e]];
      if (j === -1 || i >= j) {
        continue;
      }
      if (color[i] === 0) {
        dinic.addEdge(i, j, inf);
      } else {
        dinic.addEdge(j, i, inf);
      }
    }
  }

  dinic.maxFlow(source, sink);
  const reach = dinic.reachable(source);
  for (let i = 0; i < k; i += 1) {
    if ((color[i] === 0 && reach[i] === 1) || (color[i] === 1 && reach[i] === 0)) {
      chosenGlobal[vertices[i]] = 1;
    }
  }

  return totalWeight;
}

function buildSmallComponent(instance, vertices, localIndex) {
  const { weights, start, adj } = instance;
  const k = vertices.length;
  setLocalIndex(vertices, localIndex);

  const localWeights = new Float64Array(k);
  const compat = Array.from({ length: k }, () => new Uint8Array(k));
  const originalAdj = Array.from({ length: k }, () => []);

  for (let i = 0; i < k; i += 1) {
    localWeights[i] = weights[vertices[i]];
  }

  for (let i = 0; i < k; i += 1) {
    const gv = vertices[i];
    for (let e = start[gv]; e < start[gv + 1]; e += 1) {
      const j = localIndex[adj[e]];
      if (j !== -1) {
        originalAdj[i].push(j);
      }
    }
  }

  const adjMatrix = Array.from({ length: k }, () => new Uint8Array(k));
  for (let i = 0; i < k; i += 1) {
    for (const j of originalAdj[i]) {
      adjMatrix[i][j] = 1;
    }
  }

  for (let i = 0; i < k; i += 1) {
    for (let j = 0; j < k; j += 1) {
      if (i !== j && adjMatrix[i][j] === 0) {
        compat[i][j] = 1;
      }
    }
  }

  clearLocalIndex(vertices, localIndex);
  return { localWeights, compat };
}

function solveSmallExact(localWeights, compat) {
  const n = localWeights.length;
  let bestWeight = 0;
  let bestSet = [];
  const current = [];

  function colorSort(cands) {
    const colorVerts = [];
    const colorMax = [];
    const ordered = [];
    const bounds = [];

    for (const v of cands) {
      let c = 0;
      for (; c < colorVerts.length; c += 1) {
        let bad = false;
        for (const u of colorVerts[c]) {
          if (compat[v][u] === 1) {
            bad = true;
            break;
          }
        }
        if (!bad) {
          break;
        }
      }
      if (c === colorVerts.length) {
        colorVerts.push([]);
        colorMax.push(0);
      }
      colorVerts[c].push(v);
      colorMax[c] = Math.max(colorMax[c], localWeights[v]);
    }

    let prefix = 0;
    for (let c = 0; c < colorVerts.length; c += 1) {
      colorVerts[c].sort((a, b) => localWeights[a] - localWeights[b]);
      prefix += colorMax[c];
      for (const v of colorVerts[c]) {
        ordered.push(v);
        bounds.push(prefix);
      }
    }
    return { ordered, bounds };
  }

  function expand(cands, currentWeight) {
    if (cands.length === 0) {
      if (currentWeight > bestWeight) {
        bestWeight = currentWeight;
        bestSet = current.slice();
      }
      return;
    }

    cands.sort((a, b) => localWeights[b] - localWeights[a]);
    const { ordered, bounds } = colorSort(cands);

    for (let idx = ordered.length - 1; idx >= 0; idx -= 1) {
      if (currentWeight + bounds[idx] <= bestWeight) {
        return;
      }
      const v = ordered[idx];
      current.push(v);
      const next = [];
      for (let j = 0; j < idx; j += 1) {
        const u = ordered[j];
        if (compat[v][u] === 1) {
          next.push(u);
        }
      }
      const newWeight = currentWeight + localWeights[v];
      if (next.length === 0) {
        if (newWeight > bestWeight) {
          bestWeight = newWeight;
          bestSet = current.slice();
        }
      } else {
        expand(next, newWeight);
      }
      current.pop();
    }
  }

  const all = [];
  for (let i = 0; i < n; i += 1) {
    all.push(i);
  }
  expand(all, 0);

  const chosen = new Uint8Array(n);
  for (const v of bestSet) {
    chosen[v] = 1;
  }
  return chosen;
}

function makeRng(seed) {
  let state = BigInt.asUintN(64, BigInt(seed));
  return function next() {
    state = BigInt.asUintN(64, state + 0x9e3779b97f4a7c15n);
    let z = state;
    z = BigInt.asUintN(64, (z ^ (z >> 30n)) * 0xbf58476d1ce4e5b9n);
    z = BigInt.asUintN(64, (z ^ (z >> 27n)) * 0x94d049bb133111ebn);
    z ^= z >> 31n;
    return Number(z & ((1n << 53n) - 1n)) / Number(1n << 53n);
  };
}

function applyHeuristicSolution(instance, vertices, localIndex, chosenGlobal, seedBase) {
  const { weights, start, adj } = instance;
  const k = vertices.length;
  setLocalIndex(vertices, localIndex);

  const localWeights = new Float64Array(k);
  const localDegrees = new Int32Array(k);
  for (let i = 0; i < k; i += 1) {
    const gv = vertices[i];
    localWeights[i] = weights[gv];
    localDegrees[i] = start[gv + 1] - start[gv];
  }

  const alphas = [0.0, 0.35, 0.7, 1.0, 1.3];
  const noises = [0.0, 0.05, 0.1];
  const runs = Math.min(10, alphas.length * noises.length);

  let bestWeight = -1;
  let bestChosen = new Uint8Array(k);

  function addVertex(state, li) {
    state.inSolution[li] = 1;
    state.totalWeight += localWeights[li];
    const gv = vertices[li];
    for (let e = start[gv]; e < start[gv + 1]; e += 1) {
      const lj = localIndex[adj[e]];
      if (lj !== -1) {
        state.blockedCount[lj] += 1;
        state.blockedWeight[lj] += localWeights[li];
      }
    }
  }

  function removeVertex(state, li) {
    state.inSolution[li] = 0;
    state.totalWeight -= localWeights[li];
    const gv = vertices[li];
    for (let e = start[gv]; e < start[gv + 1]; e += 1) {
      const lj = localIndex[adj[e]];
      if (lj !== -1) {
        state.blockedCount[lj] -= 1;
        state.blockedWeight[lj] -= localWeights[li];
      }
    }
  }

  function greedyConstruct(alpha, noise, rng) {
    const order = Array.from({ length: k }, (_, idx) => idx);
    const score = new Float64Array(k);
    for (let i = 0; i < k; i += 1) {
      const perturb = 1.0 + noise * (rng() - 0.5);
      score[i] = (localWeights[i] * perturb) / Math.pow(localDegrees[i] + 1, alpha);
    }
    order.sort((a, b) => {
      if (score[b] !== score[a]) {
        return score[b] - score[a];
      }
      return localWeights[b] - localWeights[a];
    });

    const state = {
      inSolution: new Uint8Array(k),
      blockedCount: new Int32Array(k),
      blockedWeight: new Float64Array(k),
      totalWeight: 0
    };

    for (const li of order) {
      if (state.inSolution[li] === 0 && state.blockedCount[li] === 0) {
        addVertex(state, li);
      }
    }
    return { state, order };
  }

  function greedyAugment(state, order) {
    for (const li of order) {
      if (state.inSolution[li] === 0 && state.blockedCount[li] === 0) {
        addVertex(state, li);
      }
    }
  }

  function heavyInsert(state) {
    let changed = false;
    const candidates = [];
    for (let i = 0; i < k; i += 1) {
      if (state.inSolution[i] === 0 && state.blockedCount[i] > 0 && state.blockedWeight[i] < localWeights[i]) {
        candidates.push(i);
      }
    }
    candidates.sort((a, b) => (localWeights[b] - state.blockedWeight[b]) - (localWeights[a] - state.blockedWeight[a]));

    for (const li of candidates) {
      if (state.inSolution[li] === 1 || state.blockedCount[li] === 0 || state.blockedWeight[li] >= localWeights[li]) {
        continue;
      }
      const neighbors = [];
      const gv = vertices[li];
      for (let e = start[gv]; e < start[gv + 1]; e += 1) {
        const lj = localIndex[adj[e]];
        if (lj !== -1 && state.inSolution[lj] === 1) {
          neighbors.push(lj);
        }
      }
      let removedWeight = 0;
      for (const lj of neighbors) {
        removedWeight += localWeights[lj];
      }
      if (removedWeight >= localWeights[li]) {
        continue;
      }
      for (const lj of neighbors) {
        removeVertex(state, lj);
      }
      if (state.blockedCount[li] === 0) {
        addVertex(state, li);
        changed = true;
      }
    }
    return changed;
  }

  for (let run = 0; run < runs; run += 1) {
    const alpha = alphas[run % alphas.length];
    const noise = noises[Math.floor(run / alphas.length) % noises.length];
    const rng = makeRng(seedBase + run * 104729);
    const { state, order } = greedyConstruct(alpha, noise, rng);
    greedyAugment(state, order);
    for (let iter = 0; iter < 6; iter += 1) {
      const changed = heavyInsert(state);
      greedyAugment(state, order);
      if (!changed) {
        break;
      }
    }
    if (state.totalWeight > bestWeight) {
      bestWeight = state.totalWeight;
      bestChosen = state.inSolution.slice();
    }
  }

  for (let i = 0; i < k; i += 1) {
    if (bestChosen[i] === 1) {
      chosenGlobal[vertices[i]] = 1;
    }
  }

  clearLocalIndex(vertices, localIndex);
}

function solveInstance(instance) {
  const { n, deg, start, adj, weights } = instance;
  const chosen = new Uint8Array(n);
  const visited = new Uint8Array(n);
  const queue = new Int32Array(n);
  const localIndex = new Int32Array(n);
  localIndex.fill(-1);

  const SMALL_EXACT_LIMIT = 60;

  for (let root = 0; root < n; root += 1) {
    if (visited[root] === 1) {
      continue;
    }

    let head = 0;
    let tail = 0;
    queue[tail++] = root;
    visited[root] = 1;
    const vertices = [];
    let degreeSum = 0;

    while (head < tail) {
      const v = queue[head++];
      vertices.push(v);
      degreeSum += deg[v];
      for (let e = start[v]; e < start[v + 1]; e += 1) {
        const to = adj[e];
        if (visited[to] === 0) {
          visited[to] = 1;
          queue[tail++] = to;
        }
      }
    }

    const edgeCount = degreeSum >> 1;
    if (edgeCount === 0) {
      for (const v of vertices) {
        chosen[v] = 1;
      }
      continue;
    }

    if (edgeCount === vertices.length - 1) {
      applyTreeSolution(instance, vertices, localIndex, chosen);
      continue;
    }

    if (vertices.length <= SMALL_EXACT_LIMIT) {
      const small = buildSmallComponent(instance, vertices, localIndex);
      const chosenLocal = solveSmallExact(small.localWeights, small.compat);
      for (let i = 0; i < vertices.length; i += 1) {
        if (chosenLocal[i] === 1) {
          chosen[vertices[i]] = 1;
        }
      }
      continue;
    }

    const color = bipartiteColoring(instance, vertices, localIndex);
    if (color !== null) {
      applyBipartiteSolution(instance, vertices, localIndex, color, chosen);
      clearLocalIndex(vertices, localIndex);
      continue;
    }

    clearLocalIndex(vertices, localIndex);
    applyHeuristicSolution(instance, vertices, localIndex, chosen, 0xC0D3D00D + root);
  }

  const picked = [];
  let totalWeight = 0;
  for (let i = 0; i < n; i += 1) {
    if (chosen[i] === 1) {
      picked.push(i + 1);
      totalWeight += weights[i];
    }
  }

  return { totalWeight, picked, chosen };
}

function validateSolution(instance, chosen) {
  const { n, weights, edgesU, edgesV, m } = instance;
  if (chosen.length !== n) {
    throw new Error("Chosen array has wrong length");
  }
  let totalWeight = 0;
  for (let i = 0; i < n; i += 1) {
    if (chosen[i]) {
      totalWeight += weights[i];
    }
  }
  for (let i = 0; i < m; i += 1) {
    if (chosen[edgesU[i]] && chosen[edgesV[i]]) {
      throw new Error(`Invalid solution: edge (${edgesU[i] + 1}, ${edgesV[i] + 1}) is inside the set`);
    }
  }
  return totalWeight;
}

function solveText(text) {
  const instance = parseInput(text);
  const result = solveInstance(instance);
  const verified = validateSolution(instance, result.chosen);
  if (Math.abs(verified - result.totalWeight) > 0.5) {
    throw new Error("Reported weight does not match verified weight");
  }
  return `${Math.round(result.totalWeight)}\n${result.picked.join(" ")}\n`;
}

if (require.main === module) {
  const input = fs.readFileSync(0, "utf8");
  process.stdout.write(solveText(input));
}

module.exports = {
  parseInput,
  solveInstance,
  validateSolution,
  solveText
};
