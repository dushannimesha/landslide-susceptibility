# Landslide labels (NASA GLC)

- Events in bbox: 58; inside grid: 55; rain-triggered: 56
- Used for labels (inside grid, rain, accuracy <= 1.0 km): **33**
- Event years: 2007–2017  |  events in forcing period (2020): 0
- Positive cells (<= 500 m): **105**, excluded buffer: 576, negatives: 16288

Events per year:

|   date |   events |
|-------:|---------:|
|   2007 |        2 |
|   2008 |        6 |
|   2009 |        5 |
|   2010 |        2 |
|   2011 |        6 |
|   2012 |        4 |
|   2013 |        9 |
|   2014 |        8 |
|   2015 |        7 |
|   2016 |        5 |
|   2017 |        4 |

> **Note:** no catalogued events in the forcing year. Labels are spatial (historical inventory 2007–2017). For event-date (temporal) validation, re-export the GEE dataset for 2007–2017 with the same columns; the pipeline accepts any year range unchanged.