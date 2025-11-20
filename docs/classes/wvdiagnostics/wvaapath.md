---
layout: default
title: wvaapath
parent: WVDiagnostics
grand_parent: Classes
nav_order: 142
mathjax: true
---

#  wvaapath

path to a WVTransform with explicit anti-aliasing


---

## Discussion
The WVTransformBoussinesq can take a very long time to
  initialize, so if explicity anti-aliasing is requested, we cache
  a copy of the transform.
