# Strict HQCC provenance of BJ A5 seeds

Strict lattice size: **962** (excludes classical contaminants 5, 16).
20 retained if present via model `towers_20` / seeds.

16 in lattice: True; 20 in lattice: True; 5 in lattice: True

## Seeds from broad search (annotated)

| a | b | both strict? | classical 20/16? | a provenance | b provenance |
|--:|--:|:---:|:---:|---|---|
| 20 | -16 | True | True | generator, = 1+(19), = -1+(21) | generator, = 1+(-17), = -1+(-15) |
| 20 | 16 | True | True | generator, = 1+(19), = -1+(21) | generator, = 1+(15), = -1+(17) |
| -55 | -88 | True | False | generator, = 1+(-56), = -1+(-54) | generator, = 1+(-89), = -1+(-87) |
| -55 | 88 | True | False | generator, = 1+(-56), = -1+(-54) | generator, = 1+(87), = -1+(89) |
| 95 | -76 | True | False | generator, = 1+(94), = -1+(96) | generator, = 1+(-77), = -1+(-75) |
| 95 | 76 | True | False | generator, = 1+(94), = -1+(96) | generator, = 1+(75), = -1+(77) |
| 95 | -532 | True | False | generator, = 1+(94), = -1+(96) | generator, = -1+(-531), = 2*(-266) |
| 95 | 532 | True | False | generator, = 1+(94), = -1+(96) | generator, = 1+(531), = --1+(531) |
| -100 | -400 | True | False | generator, = 1+(-101), = -1+(-99) | generator, = -1+(-399), = 2*(-200) |
| -100 | 400 | True | False | generator, = 1+(-101), = -1+(-99) | generator, = 1+(399), = --1+(399) |
| 124 | -496 | True | False | generator, = 1+(123), = -1+(125) | generator, = -1+(-495), = 2*(-248) |
| 124 | 496 | True | False | generator, = 1+(123), = -1+(125) | generator, = 1+(495), = --1+(495) |
| -180 | -432 | False | False | generator, = 2*(-90), = 2+(-182) | = 2*(-216), = 2+(-434), = -2+(-430) |
| -180 | 432 | False | False | generator, = 2*(-90), = 2+(-182) | = 2*(216), = 2+(430), = -2+(434) |
| 220 | -528 | False | False | generator, = 2*(110), = -2+(222) | = 2+(-530), = --2+(-530), = 3+(-531) |
| 220 | 528 | False | False | generator, = 2*(110), = -2+(222) | = -2+(530), = -3+(531), = -4+(532) |

## Fresh strict-lattice search
- tested pairs: 272484
- even seeds: 22
- A5 seeds: **14**
- classical (20,±16) among them: 2

### Non-classical strict A5 seeds

- a=-55 b=-88: `x**5 - 55*x - 88` Gal=S5TransitiveSubgroups.A5 prov_a=['generator', '= 1+(-56)', '= -1+(-54)', '= --1+(-56)'] prov_b=['generator', '= 1+(-89)', '= -1+(-87)', '= --1+(-89)']
- a=-55 b=88: `x**5 - 55*x + 88` Gal=S5TransitiveSubgroups.A5 prov_a=['generator', '= 1+(-56)', '= -1+(-54)', '= --1+(-56)'] prov_b=['generator', '= 1+(87)', '= -1+(89)', '= --1+(87)']
- a=95 b=-76: `x**5 + 95*x - 76` Gal=S5TransitiveSubgroups.A5 prov_a=['generator', '= 1+(94)', '= -1+(96)', '= --1+(94)'] prov_b=['generator', '= 1+(-77)', '= -1+(-75)', '= --1+(-77)']
- a=95 b=76: `x**5 + 95*x + 76` Gal=S5TransitiveSubgroups.A5 prov_a=['generator', '= 1+(94)', '= -1+(96)', '= --1+(94)'] prov_b=['generator', '= 1+(75)', '= -1+(77)', '= --1+(75)']
- a=95 b=-532: `x**5 + 95*x - 532` Gal=S5TransitiveSubgroups.A5 prov_a=['generator', '= 1+(94)', '= -1+(96)', '= --1+(94)'] prov_b=['generator', '= -1+(-531)', '= 2*(-266)', '= -2+(-530)']
- a=95 b=532: `x**5 + 95*x + 532` Gal=S5TransitiveSubgroups.A5 prov_a=['generator', '= 1+(94)', '= -1+(96)', '= --1+(94)'] prov_b=['generator', '= 1+(531)', '= --1+(531)', '= 2*(266)']
- a=-100 b=-400: `x**5 - 100*x - 400` Gal=S5TransitiveSubgroups.A5 prov_a=['generator', '= 1+(-101)', '= -1+(-99)', '= --1+(-101)'] prov_b=['generator', '= -1+(-399)', '= 2*(-200)', '= -2*(200)']
- a=-100 b=400: `x**5 - 100*x + 400` Gal=S5TransitiveSubgroups.A5 prov_a=['generator', '= 1+(-101)', '= -1+(-99)', '= --1+(-101)'] prov_b=['generator', '= 1+(399)', '= --1+(399)', '= 2*(200)']
- a=124 b=-496: `x**5 + 124*x - 496` Gal=S5TransitiveSubgroups.A5 prov_a=['generator', '= 1+(123)', '= -1+(125)', '= --1+(123)'] prov_b=['generator', '= -1+(-495)', '= 2*(-248)', '= -2+(-494)']
- a=124 b=496: `x**5 + 124*x + 496` Gal=S5TransitiveSubgroups.A5 prov_a=['generator', '= 1+(123)', '= -1+(125)', '= --1+(123)'] prov_b=['generator', '= 1+(495)', '= --1+(495)', '= 2*(248)']
- a=320 b=-512: `x**5 + 320*x - 512` Gal=S5TransitiveSubgroups.A5 prov_a=['generator', '= 2*(160)', '= 2+(318)', '= -2+(322)'] prov_b=['generator', '= 1+(-513)', '= --1+(-513)', '= 2*(-256)']
- a=320 b=512: `x**5 + 320*x + 512` Gal=S5TransitiveSubgroups.A5 prov_a=['generator', '= 2*(160)', '= 2+(318)', '= -2+(322)'] prov_b=['generator', '= -1+(513)', '= 2*(256)', '= -2*(-256)']

### Homogenised non-classical strict families

#### seed (-55, -88)
- proved_even=True n_A5=10
- family `x**5 + (-55)*t**4*x + (-88)*t**5`
- groups `{'S5TransitiveSubgroups.A5': 10}`

#### seed (-55, 88)
- proved_even=True n_A5=10
- family `x**5 + (-55)*t**4*x + (88)*t**5`
- groups `{'S5TransitiveSubgroups.A5': 10}`

#### seed (95, -76)
- proved_even=True n_A5=10
- family `x**5 + (95)*t**4*x + (-76)*t**5`
- groups `{'S5TransitiveSubgroups.A5': 10}`

#### seed (95, 76)
- proved_even=True n_A5=10
- family `x**5 + (95)*t**4*x + (76)*t**5`
- groups `{'S5TransitiveSubgroups.A5': 10}`

#### seed (95, -532)
- proved_even=True n_A5=10
- family `x**5 + (95)*t**4*x + (-532)*t**5`
- groups `{'S5TransitiveSubgroups.A5': 10}`

#### seed (95, 532)
- proved_even=True n_A5=10
- family `x**5 + (95)*t**4*x + (532)*t**5`
- groups `{'S5TransitiveSubgroups.A5': 10}`

#### seed (-100, -400)
- proved_even=True n_A5=10
- family `x**5 + (-100)*t**4*x + (-400)*t**5`
- groups `{'S5TransitiveSubgroups.A5': 10}`

#### seed (-100, 400)
- proved_even=True n_A5=10
- family `x**5 + (-100)*t**4*x + (400)*t**5`
- groups `{'S5TransitiveSubgroups.A5': 10}`

#### seed (124, -496)
- proved_even=True n_A5=10
- family `x**5 + (124)*t**4*x + (-496)*t**5`
- groups `{'S5TransitiveSubgroups.A5': 10}`

#### seed (124, 496)
- proved_even=True n_A5=10
- family `x**5 + (124)*t**4*x + (496)*t**5`
- groups `{'S5TransitiveSubgroups.A5': 10}`

#### seed (320, -512)
- proved_even=True n_A5=10
- family `x**5 + (320)*t**4*x + (-512)*t**5`
- groups `{'S5TransitiveSubgroups.A5': 10}`

#### seed (320, 512)
- proved_even=True n_A5=10
- family `x**5 + (320)*t**4*x + (512)*t**5`
- groups `{'S5TransitiveSubgroups.A5': 10}`

## Verdict

Found **12** strict-HQCC A5 BJ seeds distinct from classical (20,±16). These admit the same homogenised theorem: disc(f_t)=t^20 disc(seed) square.

_Generated by hqcc_strict_analysis.py_