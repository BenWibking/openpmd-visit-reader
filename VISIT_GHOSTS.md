# Supporting VisIt Ghost Synthesis

VisIt can synthesize ghost zones for structured or AMR meshes when a database
plugin delivers three pieces of auxiliary information:

1. **Structured domain boundaries** – an `avtStructuredDomainBoundaries`
   instance populated with logical indices and levels for each patch.
2. **Global node IDs** – deterministic IDs for every mesh node in a patch,
   returned as a `vtkIdTypeArray` named `avtGlobalNodeId` when the engine asks
   for `GLOBAL_NODE_IDS` / `AUXILIARY_DATA_GLOBAL_NODE_IDS`.
3. **Global zone IDs** – similar IDs for every cell, returned as
   `avtGlobalZoneId` in response to `GLOBAL_ZONE_IDS` queries.

Without this information VisIt assumes there are no overlapping zones between
patches, so the `InverseGhostZone` operator ends up stripping all data. The
`avtGhostZoneAndFacelistFilter` also explicitly demands both global-ID arrays
before it will trust domain adjacency.

### Implementation sketch

* Cache an `avtStructuredDomainBoundaries` per mesh (usually while populating
  metadata) via `cache->CacheVoidRef(meshName, AUXILIARY_DATA_DOMAIN_BOUNDARY_INFORMATION, timestep, -1, vr)`. VisIt never explicitly asks for `domain == -1` later; it only retrieves what the plugin tucked away in the cache. Serving local `avtLocalStructuredDomainBoundaryList`s for individual domains not only doesn’t help, it can make VisIt reject the entire boundary set (especially in 2‑D where the unused axis still contributes nodes). Once the global object is cached, simply let the base `avtMTMDFileFormat::GetAuxiliaryData` return it. This is exactly what the built-in AMR readers (SAMRAI, FLASH, BATL, Chombo, etc.) do.

* Build deterministic global IDs by flattening logical indices into a single
  `vtkIdType`. For a cell patch the formula is
  `gi + Nx * (gj + Ny * gk)` where `gi/gj/gk` are logical coordinates and
  `Nx/Ny` are the mesh-wide dimensions along x/y. Node-centered meshes skip the
  extra `+1` on the upper bound; zone-centered meshes add it. Returning freshly
  allocated arrays works, but caching them in VisIt’s `VoidRefCache` (as the
  SAMRAI reader does) avoids recomputation on every request and guarantees the
  lifetime matches VisIt’s expectations.

  *Each AMR level shares this same global logical index space.* Patch
  `logicalLower/logicalUpper` values are already expressed on the finest grid,
  so finer levels simply produce higher index ranges; there’s no practical
  limit on the number of AMR levels beyond the 32-bit logical-extent ceiling
  VisIt imposes for structured meshes.

* Continue advertising `AVT_HAS_GHOSTS` in mesh metadata, but do **not** try to
  fabricate `avtGhostZones` arrays yourself. Once VisIt accepts the structured
  boundaries and global IDs it will mark overlap regions as ghosts, allowing
  operators like `InverseGhostZone` to work.

### Debugging tips

* Run VisIt with `-debug 5` to capture `*.vlog` logs in the repo root. Look for
  `[openpmd-api-plugin] Global node ids ready...` to confirm the plugin served
  the arrays, and for `Rejecting domain boundaries...` messages if VisIt still
  distrusts the structured object.

* The rejection warning also prints the point count VisIt expected versus what
  the dataset actually has. If they differ, double-check the extents passed to
  `SetIndicesForAMRPatch`—they must be inclusive-exclusive (lower, upper+1).

* When things work, the logs show `avtGenericDatabase::GetDomainBoundaryInformation: Global for local success` with no subsequent rejection, and the
  `InverseGhostZone` operator stops discarding the entire mesh.
