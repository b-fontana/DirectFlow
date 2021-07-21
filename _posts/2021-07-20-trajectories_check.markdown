---
layout: post
title:  "Consistency check: x=5cm, y=5cm"
date:   2021-07-20 10:00:18 +0200
categories: plots
---

We deactivate the two magnetic dipoles along the ```X``` direction (the ones which introduced an assymetry between positive and negative ```Z```, labelled "D_corr" and "Muon"). We observe that the "kick" along the ```Y``` axis is thus gone. In addition, by starting both beams at exactly mirrored positions, we show that they intersect in the middle (```Z```=0), as expected.

-------------

{% include traj_libs.html %}

{% include traj_5_5_1380.html %}

-----------

You can find the code used for plotting [here][plotcode].

[plotcode]: https://github.com/b-fontana/DirectFlow/blob/master/python/trajectory.py
