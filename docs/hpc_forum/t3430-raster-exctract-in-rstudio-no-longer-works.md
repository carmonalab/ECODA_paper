# Raster::exctract in Rstudio no longer works

- Source: https://hpc-community.unige.ch/t/3430

- Created: 2024-04-30T07:33:12.160Z

- Posts: 4

- Category: 14

- Pinned: False

- Snapshot: 2026-08-11

---

## Post 1 by @Mariana.Bouvier-Rabe (2024-04-30T07:33:12.214Z)

Hello,
last week I ran this command:
```
raster_altitude_agreg ← aggregate(edm, fact=res(tasr1)/res(edm), fun=mean, na.rm=TRUE)
altitudeMoyenne ← raster::extract(raster_altitude_agreg , v_elev)
```
And it worked perfectly.
But now, it doesn t work anymore. I still get an error message:
```
AltitudeMoyenne ← raster::extract(raster_altitude_agreg,coord_lonlatC)

Error in file(fn, “rb”) : cannot open the connection
In addition: Warning message:
In file(fn, “rb”) :
cannot open file ‘/tmp/RtmpiwjSJ7/raster/r_tmp_2024-04-29_084403.028646_2725427_25595.gri’: No such file or directory

altitudeMoyenne ← raster::extract(raster_altitude_agreg , v_elev)
Error in file(fn, “rb”) : cannot open the connection
In addition: Warning message:
In file(fn, “rb”) :
cannot open file ‘/tmp/RtmpiwjSJ7/raster/r_tmp_2024-04-29_084403.028646_2725427_25595.gri’: No such file or directory
```
I am not sure to understand what is happening ?
My working directory is normally the same.
This is the same for writeRaster() function.
Thank you in advance.
Mariana


## Post 2 by @Yann.Sagon (2024-05-02T13:38:51.147Z)

Dear @Mariana.Bouvier-Rabe[@Mariana.Bouvier-Rabe](https://hpc-community.unige.ch/u/mariana.bouvier-rabe)  the error message seems related with your code, not Rstudio. Try to share your code with us, maybe someone with R knowledge will be able to help you.
Best


## Post 3 by @Mariana.Bouvier-Rabe (2024-05-03T17:50:19.205Z)

Hello,
Thank you for your answer. Finally, the problem no longer occurred. I did’nt understand what happened.
Have a nice day,
Mariana


## Post 4 by @Matthieu.Stigler (2024-06-06T08:35:11.671Z)

In case this happens again: this could be due to the fact that  `raster` is using temporary files when the raster doesn’t fit in memory. So you would need to look at functions like  `showTmpFiles()`, eventually setting your own `rasterOptions(tmpdir = "D:/rtmp/")`, to try to understand where the problem is.
Note also that `raster` has now been superseded by `terra`, so ideally use `terra::extract()` instead?
