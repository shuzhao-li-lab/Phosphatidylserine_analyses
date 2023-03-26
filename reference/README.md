# reference

HI Minghao,
 
I have completed the analysis generating level-4 annotations using the computationally generated PS class lipids. I was able to generate PS, LPS, LPS O-, and LPS O- ;O variants of PS lipids. After deduplicating isomeric formulas, this yields 3129 unique formulas. I searched only for the [M-H+e] adducts and I used the mpmath library in python to ensure the neutral masses and adduct masses were accurate.
 
I incorporated the improvements you suggested during our meeting today. Now the annotations are generated in a two-step process. First the m/z of the monoisotopic (e.g., m+13C0) isotopologue is searched against the feature list with a 5ppm tolerance. In the second step, the 13C isotopologues of the formula are generated in ascending order and searched for in the feature list by m/z with a 5ppm tolerance and with a requirement that the retention time of possible isotopologues be within the rtime_right_base and rtime_left_base of the m+13C0 feature with some wiggle room, rtime - rtime_left_base on the left and ritme_right_base – time on the right to allow for some rtime error. This search continues until no such isotopolgue is found. Since the number of isotopologues observed is highly dependent on feature intensity, simply having more isotopolgues observed does not imply the assignment is necessarily correct but it is some circumstantial evidence that it is. This is recorded in the “hits_isotopologue_chain” field of the output json. I did not use any intensity information.
 
Some lipids have up to 4 13C isotopologues identified which is subjective evidence that this is a plausible assignment.
 
I’ve double, triple checked the formulas my code generates against selected entries in LipidMaps and on the target list you had, and they appear correct. I would appreciate you double checking any possible hits you report to ensure I did not make a miscalculation.
 
I need to get this code on GitHub so we can reuse it in the future, but I want to figure out the GitHub organization nonsense first (hopefully tomorrow).
 
Here is an example output with descriptions:
 
        "name": "PS 38:0", # common name
        "neutral_mass": 819.59893497033, # neutral mass
        "formula": "C44H86NO10P", #formula
        "formula_dict": { # dict representation of formula
            "C": 44,
            "H": 86,
            "O": 10,
            "N": 1,
            "P": 1
        },
        "isomers": [], # common names of isomers for this lipid (if there any)
        "[M-H+e]": 818.591658518009, # mass of M-H+e adduct for searching
        "mz_only_hits": [ # hits from the first round of searching, m/z only
            "F12441",
            "F12442"
        ],
        "hits_isotopologue_chain": { # isotopologues found for each mz only hit
            "F12442": {
                "F12442": {
                    "0": "F12442" #m+13C0 for F12442
                },
                "1": [
                    "F12548" #m+13C1 for F12442
                ],
                "2": [
                    "F12676" #m+13C2 for F12442
                ],
                "3": [
                    "F12787" #m+13C3 for F12442
                ],
                "4": [
                    "F12921" #m+13C4 for F12442
                ]
            }
        }
    },
 