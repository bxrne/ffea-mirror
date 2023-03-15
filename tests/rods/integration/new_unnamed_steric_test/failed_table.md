| Test name                | Config OK | FFEA finish | stderr |  Last step  | Note |
|--------------------------|-----------|-------------|-------|-------------|-------       |
| oblique_Cross+x+RxSmall  |    y       |    n     |   Sum of node forces not equal to element force.    |   0    |       |
| oblique_Cross+x-RxSmall  |    y       |    n     |   Sum of node forces not equal to element force.    |   0    |       |
| oblique_Cross+y+RySmall  |     y      |    n     |   Sum of node forces not equal to element force.    |   0    |       |
| oblique_Cross+y-RySmall  |     y      |    n     |   Sum of node forces not equal to element force.    |   0    |       |
| oblique_Cross-x+RxSmall  |     y      |    n     |    Sum of node forces not equal to element force.   |   0   |       |
| oblique_Cross-x-RxSmall  |     y      |    n     |    Sum of node forces not equal to element force.   |   0    |       |
| oblique_Cross-y+RySmall  |     y      |    n     |    Sum of node forces not equal to element force.   |   0    |       |
| oblique_Cross-y-RySmall  |     y      |    n     |    Sum of node forces not equal to element force.   |   0    |       |
| oblique_T+z+RxBig        |     y      |    n     |  Sum of node forces not equal to element force. |   0    |       |
| oblique_zPlane+y+RxSmall |     n      |    n     |  Sum of node forces not equal to element force. |   0    |    Rods should intersect at end, but instead do in the middle, and fully through   |
| oblique_zPlane+y+RxBig |     n      |    n     |  Sum of node forces not equal to element force. |   0    |    See above cell.  |
| parallel_+xyz            |    y       |    y     |   n    |   5000    |  no steric force vectors visible
| parallel_+xzHalf         |     y      |    n     |  Sum of node forces not equal to element force. |   0    |       |
| parallel_+z              |      y     |    y     |   n    |   5000    | no steric force vectors visible
| parallel_-xyz            |      y     |    y     |   n    |   5000    | no steric force vectors visible
| parallel_-xzHalf         |     y      |    n     |  Sum of node forces not equal to element force. |   0    |       |
| parallel_-z              |      y     |    y     |   n    |   5000    | no steric force vectors visible
| perp_L+x+z+Ry            |      y     |    y     |   n    |   5000    |  repulsion occurs
| perp_L+y+z+Rx            |     y      |    y     |   n    |   5000    |  repulsion occurs
| perp_L+y-z+Rx            |     n*      |    y     |   n    |   5000    |  incorrect entry in config creation script. *now fixed
| perp_L-x+z+Ry            |     y      |    y     |   n    |   5000    |  repulsion occurs
| perp_L-x-z+Ry            |      y     |    n     |  Sum of node forces not equal to element force. |    3300   |       |
| perp_L-y+z+Rx            |     y      |    y     |   n    |   5000    | repulsion occurs
| perp_L-y-z+Rx            |     n*      |    y     |   n    |   5000    | repulsion occurs. *now fixed
