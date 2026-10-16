# Source provenance and license status

This local branch combines the supplied KSSOLV tree and the earlier WI source release. It does not claim sole WI authorship of KSSOLV, solver dependencies, upstream Quantum ESPRESSO routines, or pseudopotentials. Original source headers and included third-party license files remain intact.

## Existing notices preserved

| Component | Evidence retained in this branch |
| --- | --- |
| KSSOLV startup and many class/solver files | Upstream MIT notices in source, including Yingzhou Li and Chao Yang. This does not establish a blanket license over every bundled file. |
| export_fig | [Original license](external/export_fig/LICENSE), copyright 2014 Oliver J. Woodford and Yair M. Altman; BSD-style three-condition notice preserved. Java source/class assets also retained. |
| HANSO | [Original GPLv3 text](external/hanso2_2/GPLver3.license) and existing source comments. |
| NLCG | [Original GPLv3 text](external/nlcg1_0/GPLver3.license) and existing source comments. |
| Screened-exchange Fortran/MEX | `src/SCF_DCM/Potential/Vxc/VGCExchange/pbexsrF.F90` retains its 2002-2009 Quantum ESPRESSO GNU GPL header. The supplied Linux binary was excluded because it embeds a local build path; source is retained for local recompilation. Exact upstream revision, GPL version/options and corresponding distribution terms need review. |
| UPF/PSP data | Supplied `ppdata` files retained with original embedded metadata. Individual provider/version/redistribution terms must be established before public distribution; filenames/hashes are not permission. |
| MATLAB/toolboxes and Wannier90 executable | External, not bundled. Their actual distributions' terms apply. |

The base has no authoritative repository-wide license or Git metadata identifying its upstream revision. No new blanket MIT/BSD/GPL license was assigned by this integration, and no `LICENSE` file was fabricated. Full provided GPLv3 texts were retained unchanged; their presence alone does not settle the WI translation's applicable GPL version.

## WI translation provenance remains unresolved

The original WI project's `README_SCDM.md`, line 6, states that MATLAB SCDM code was translated from Quantum ESPRESSO's `pw2wannier90.f90`; lines 273-275 identify GNU GPL. The original reference Fortran header credits the Quantum ESPRESSO/Wannier90 groups (2003-2013), with SCDM contributions by Valerio Vitale and noncollinear SCDM by Jae-Mo Lihm, among other named authors. A manuscript program-summary field instead lists BSD 3-clause. These are historical evidence from the source project, not a resolved license grant for this branch.

The upstream `pw2wannier90.f90` reference was not copied into the earlier WI package or this branch. Omitting that reference does not remove obligations attached to translated/derived MATLAB code. Before publication, identify copied/translated portions, their source revisions and applicable license versions, retain all required attributions/texts, and obtain an authorized compatible project licensing decision. Reconcile the manuscript field with the actual software terms.

Merging KSSOLV and WI makes the runtime source tree self-contained for the documented checks; it does not resolve the GPL/BSD issue or authorize PP redistribution. Public release remains pending these decisions. The provided screened-exchange GPL source and GPLv3 dependencies are additional reasons not to describe the whole tree as solely MIT/BSD.

## Runtime records

The unchanged WI provenance writer can snapshot solver source, PP data, adapted SOC classes and local absolute paths into case results. Such snapshots are ignored and were removed from this delivered source folder after validation. Review rights, private paths and size before later uploading generated archives/checkpoints. All sources/data in the manifest retain their component identity; they are not relabeled as original WI code.
