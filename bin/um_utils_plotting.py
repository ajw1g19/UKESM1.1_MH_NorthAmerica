#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import regionmask as rm
from shapely.geometry import Polygon, LineString, Point
from shapely.ops import split, unary_union
from shapely.geometry.polygon import orient
import xarray as xr

def plotBaseMap(ax, lims=None,
                coastres="110m", coastlw=0.3, coastcol="grey",
                borders=True, borderlw=0.15, bordercol="grey",
                land=False, landcol="lightgrey", landalph=1.0,
                ocean=False, oceancol="white",
                lakes=False, lakelw=0.3, lakecol="lightblue",
                states=False, statelw=0.1, statelinecol="grey", statescale="50m", statecol="none",
                spinelw=1, spinecol="black"):
    """
    Plots a basemap on a given geo axes with coastlines and country borders.
    Allows control over coastline width, border width, colors, ocean fill, and map extent.
    Args:
        ax: matplotlib/cartopy geo axes
        lims: map extent (list or None)
        coastlw: coastline line width
        coastcol: coastline color
        borders: whether to plot country borders
        borderlw: border line width
        bordercol: border color
        ocean: whether to fill ocean
        oceancol: ocean fill color
        spinelw: geo spine line width
        spinecol: geo spine color
    Returns:
        None
    """    
    ax.spines['geo'].set_edgecolor(spinecol)
    ax.spines["geo"].set_linewidth(spinelw)
    ax.coastlines(resolution=coastres, lw=coastlw, edgecolor=coastcol)
    
    if borders:
        ax.add_feature(cfeature.BORDERS, edgecolor=bordercol, linewidth=borderlw, zorder=2)

    if land:
        ax.add_feature(cfeature.LAND, facecolor=landcol, alpha=landalph, zorder=1)

    if ocean:
        ax.add_feature(cfeature.OCEAN, facecolor=oceancol, zorder=1)

    if lakes:
        ax.add_feature(cfeature.LAKES, edgecolor="black", facecolor=lakecol, linewidth=lakelw, zorder=2)

    if states:
        ax.add_feature(cfeature.NaturalEarthFeature(category="cultural", name="admin_1_states_provinces_lines", scale=statescale, facecolor=statecol, edgecolor=statelinecol), linewidth=statelw, zorder=1)

    try:  
        ax.set_extent(lims, crs=ccrs.PlateCarree())         
    except: 
        ax.set_global()



def plotGrid(ax, xlabeledlines, ylabeledlines, xlines=None, ylines=None, 
             labelsize=8, lw=0.4, ls=(0,(5,10)), col="k", alph=0.8, 
             top_lbl=False, right_lbl=False):
    """
    Adds labelled and non-labelled gridlines to a geo axes.
    Allows customization of gridline locations, label size, line width, style, color, alpha, and label placement.
    Args:
        ax: matplotlib/cartopy geo axes
        xlabeledlines: x gridline locations with labels
        ylabeledlines: y gridline locations with labels
        xlines: x gridline locations without labels
        ylines: y gridline locations without labels
        labelsize: label font size
        lw: gridline width
        ls: gridline linestyle
        col: gridline color
        alph: gridline alpha
        top_lbl: show top labels
        right_lbl: show right labels
    Returns:
        None
    """    
    # Labels only - always vector
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, lw=0, alpha=0.0)
    gl.top_labels, gl.right_labels = top_lbl, right_lbl
    gl.xlabel_style = {"size": labelsize}
    gl.ylabel_style = {"size": labelsize}
    gl.xlocator = mticker.FixedLocator(xlabeledlines)
    gl.ylocator = mticker.FixedLocator(ylabeledlines)
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # Lines only - rasterizable without touching any text
    glines = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                          lw=lw, color=col, alpha=alph, linestyle=ls)
    glines.xlocator = mticker.FixedLocator(list(xlabeledlines) + list(xlines or []))
    glines.ylocator = mticker.FixedLocator(list(ylabeledlines) + list(ylines or []))

    return [gl, glines]



def find_axis(fig, title=None, xlabel=None, ylabel=None, contains=False):
    """
    Find an Axes in fig by title/xlabel/ylabel.
    Matches center/left/right titles and any Text artists.
    If contains=True does substring matching.
    Returns first matching Axes or None.
    """
    def _matches(s, target):
        if s is None:
            return False
        s = str(s)
        return (target in s) if contains else (s == target)

    for ax in fig.get_axes():
        if title is not None:
            # center title
            if _matches(ax.get_title() or "", title):
                return ax

            # left/right title attributes (matplotlib stores these on the Axes)
            for attr in ("_left_title", "_right_title"):
                t = getattr(ax, attr, None)
                try:
                    txt = t.get_text() if t is not None else ""
                except Exception:
                    txt = ""
                if _matches(txt, title):
                    return ax

            # any Text artists attached to the Axes (covers manual text/title variants)
            for txt_obj in getattr(ax, "texts", []):
                try:
                    txt = txt_obj.get_text()
                except Exception:
                    continue
                if _matches(txt, title):
                    return ax

            # fallback: check all children for Text instances
            for ch in ax.get_children():
                if isinstance(ch, mpl.text.Text):
                    try:
                        txt = ch.get_text()
                    except Exception:
                        continue
                    if _matches(txt, title):
                        return ax

        if xlabel is not None:
            if _matches(ax.get_xlabel() or "", xlabel):
                return ax

        if ylabel is not None:
            if _matches(ax.get_ylabel() or "", ylabel):
                return ax

    return None


def ax_adjust(ax, shrinkx=1.0, shrinky=1.0, shiftx=0.0, shifty=0.0):
    """
    Adjust position of a single Axes or an iterable of Axes.
    Args:
        ax: matplotlib Axes instance or iterable of Axes
        shrinkx/shrinky: scale factors for width/height
        shiftx/shifty: shifts to apply to the center (in axes fraction)
    """
    def _adjust_single(a):
        pos = a.get_position()
        new_h = pos.height * shrinky
        new_w = pos.width * shrinkx

        tmp_y = pos.y0 + (pos.height - new_h) / 2.0
        tmp_x = pos.x0 + (pos.width - new_w) / 2.0

        new_x = tmp_x + shiftx
        new_y = tmp_y + shifty

        a.set_position([new_x, new_y, new_w, new_h])

    if hasattr(ax, "get_position"):
        _adjust_single(ax)
        return

    try:
        for a in ax:
            if hasattr(a, "get_position"):
                _adjust_single(a)
    except TypeError:
        raise TypeError("ax must be a matplotlib Axes or an iterable of Axes")


def preDefinedRegions():
    """
    Return lists of polygons, names, and abbreviations for pre-defined regions.
    Returns:
        raw_polys: list of polygons (list of [lon, lat] pairs)
        names: array of region names
        abbrevs: array of region abbreviations
    """
    raw_polys = [
        [[-100, 33], [-92, 33], [-92, 55], [-115, 55]], # GreatPlains
        [[-105, 25], [-90, 25], [-90, 35], [-105, 35]], # GPLLJ
        [[-168, 25], [-50, 25], [-50, 75], [-168, 75]], # N_America
    ]

    names = np.array([
        "Great Plains",
        "Great Plains Low Level Jet",
        "North America"
    ])
    
    abbrevs = np.array([
        "GP",
        "GPLLJ",
        "NAM"
    ])

    return raw_polys, names, abbrevs



def umRegion(abbrev_or_coords, lon360=True, name=None, abbrev=None):
    """
    Create a regionmask.Regions object for a specified region abbreviation or custom coordinates.
    Checks if the entered abbreviation is valid for pre-defined regions.
    Converts longitudes as needed and handles polygons crossing the prime meridian or dateline.
    Args:
        abbrev_or_coords: region abbreviation (string) or list of [lon, lat] pairs
        lon360: bool, if True output longitudes in 0:360, else -180:180
        name: region name (required for custom coordinates)
        abbrev: region abbreviation (required for custom coordinates)
    Returns:
        regionmask.Regions object
    """  
    def convert_poly(raw_poly, lon360):
        """
        Convert a polygon defined by longitude-latitude pairs to ensure all longitudes are in the 0:360 or -180:180 range.
        If the polygon crosses the prime meridian or dateline, split it appropriately.
        Adjust longitude values by adding or subtracting 360 as needed.
        Args:
            raw_poly: list of [longitude, latitude] pairs defining the polygon
            lon360: bool, if True output longitudes in 0:360, else -180:180
        Returns:
            poly_list: list of polygons with longitudes in the requested range
        """
        poly_obj = Polygon(raw_poly)
        min_lon = min([p[0] for p in raw_poly])
        max_lon = max([p[0] for p in raw_poly])

        crosses_prime_meridian = (min_lon < 0) and (max_lon > 0)
        crosses_dateline = (max_lon - min_lon) > 180

        # If input coords are 0:360 range, any polygon crossing the prime meridian needs splitting in two
        if lon360: 
            if crosses_dateline and not crosses_prime_meridian:
                unwrapped = [[lon + 360 if lon < 0 else lon, lat] for lon, lat in raw_poly]
                poly_obj = orient(Polygon(unwrapped), sign=1.0)
                split_line = LineString([(180, -90), (180, 90)])
                split_result = split(poly_obj, split_line)
                poly_list = [np.array(p.exterior.coords) for p in split_result.geoms]            
            elif crosses_prime_meridian: 
                split_line = LineString([(0, -90), (0, 90)])
                split_result = split(poly_obj, split_line)
                poly_list = [np.array(p.exterior.coords) for p in split_result.geoms]
            else:
                poly_list = [np.array(poly_obj.exterior.coords)]        
                
            for i, poly in enumerate(poly_list):
                new_poly = []
                converted_lon = any(lon < 0 for lon, lat in poly)
                for lon, lat in poly:
                    new_lon = lon + 360 if lon < 0 else lon
                    if lon == 0 and converted_lon:
                        new_lon = 360
                    new_poly.append([new_lon, lat])
                poly_list[i] = new_poly

        # If input coords are -180:180 range, any polygon crossing the dateline needs splitting in two
        elif not lon360:
            if crosses_dateline and not crosses_prime_meridian:
                unwrapped_poly = [[lon + 360 if lon < 0 else lon, lat] for lon, lat in raw_poly]
                poly_obj = orient(Polygon(unwrapped_poly), sign=1.0)
                split_line = LineString([(180, -90), (180, 90)])
                split_result = split(poly_obj, split_line)
                poly_list = [np.array(p.exterior.coords) for p in split_result.geoms]

                for i, poly in enumerate(poly_list):
                    new_poly = []
                    converted_lon = any(lon > 180 for lon, lat in poly)
                    for lon, lat in poly:
                        new_lon = lon - 360 if lon > 180 else lon
                        if lon == 180 and converted_lon:
                            new_lon = -180
                        new_poly.append([new_lon, lat])
                    poly_list[i] = new_poly

            else:
                poly_list = [np.array(poly_obj.exterior.coords)]

        return poly_list

    if isinstance(abbrev_or_coords, str):
        raw_polys, names, abbrevs = preDefinedRegions()
        if abbrev_or_coords not in abbrevs:
            raise ValueError(f"'{abbrev_or_coords}' is not a valid region abbreviation. Valid options are: {list(abbrevs)}")            
        poly_idx = np.where(abbrevs == abbrev_or_coords)[0][0]
        poly_list = convert_poly(raw_polys[poly_idx], lon360)

        n_parts = len(poly_list)
        name_list = [names[poly_idx]] * n_parts
        abbrev_list = [abbrevs[poly_idx]] * n_parts
        
        Region = rm.Regions(poly_list, names=name_list, abbrevs=abbrev_list)

    elif isinstance(abbrev_or_coords, (list, np.ndarray)):
        if name is None or abbrev is None:
            raise ValueError("For custom coordinates, you must supply both 'name' and 'abbrev'.")
        
        poly_list = convert_poly(abbrev_or_coords, lon360)
        Region = rm.Regions(poly_list, names=name, abbrevs=abbrev)

    return Region



def umRegionMask(mask, longitude, latitude, category=None, rgname="CustomRegion", rgabbrev="Custom", **kwargs):
    """
    Generate a region mask for a specified region name and longitude/latitude arrays.
    Supports user-defined, SREX, and AR6 regions, and optional land/ocean masking.
    Args:
        mask: region name (string)
        longitude: longitude array
        latitude: latitude array
        category: region category ('srex', 'ar6', or None)
        kwargs: optional 'land' argument for land/ocean masking
    Returns:
        Region_mask: boolean mask (xarray.DataArray)
    """   
    if category is None:
        if isinstance(mask, str):
            if np.max(longitude) > 180:
                Regions = umRegion(mask)
            else:
                Regions = umRegion(mask, lon360=False)

        elif isinstance(mask, (list, np.ndarray)):
                if np.max(longitude) > 180:
                    Regions = umRegion(mask, name=rgname, abbrev=rgabbrev)
                else:
                    Regions = umRegion(mask, lon360=False, name=rgname, abbrev=rgabbrev)    

        else:
            raise ValueError("mask must be a string or a list of coordinates when category is None.")
       

        masks = Regions.mask_3D(longitude, latitude)
        if len(Regions.numbers) > 1:
            indices = np.arange(len(Regions.numbers))
            Region_mask = masks.isel(region=indices).any(dim="region").squeeze()
        else:
            Region_mask = masks.squeeze()

        if "lat" and "lon" in Region_mask.dims:
            Region_mask = Region_mask.rename({"lat": "latitude", "lon": "longitude"})
    
    elif category == "srex":
        if not isinstance(mask, str):
            raise ValueError("For SREX regions, mask must be a string abbreviation.")
        valid_abbrevs = list(rm.defined_regions.srex.abbrevs)
        if mask not in valid_abbrevs:
            raise ValueError(f"'{mask}' is not a valid SREX region abbreviation. Valid options are: {valid_abbrevs}")        
        regions = rm.defined_regions.srex
        masks = regions.mask(longitude, latitude)
        Region_mask = masks.cf == mask

    elif category == "ar6":
        if not isinstance(mask, str):
            raise ValueError("For AR6 regions, mask must be a string abbreviation.")
        valid_abbrevs = list(rm.defined_regions.ar6.all.abbrevs)
        if mask not in valid_abbrevs:
            raise ValueError(f"'{mask}' is not a valid AR6 region abbreviation. Valid options are: {valid_abbrevs}")        
        regions = rm.defined_regions.ar6.all
        masks = regions.mask(longitude, latitude)
        Region_mask = masks.cf == mask
    
    if "land" in kwargs.keys():
        if "lsm" not in kwargs:
            raise ValueError("If 'land' is specified, you must also provide a 2D 'lsm' array (land=1, ocean=0).")
        lsm = kwargs["lsm"]
        lsm_bool = np.asarray(lsm) == 1

        if kwargs["land"] == True:  # Land points only
            Region_mask = np.logical_and(Region_mask, lsm_bool)
        elif kwargs["land"] == False:  # Ocean points only
            Region_mask = np.logical_and(Region_mask, ~lsm_bool)
        else:
            raise ValueError("'land' kwarg must be True or False.")

    return Region_mask



def proxyRegionMask(proxy_vals, region, lon_arr=None, lat_arr=None, landsea_mask=None):
    # Convert region polygon coordinates to 0:360 longitude range
    def _poly_to_360(poly):
        # convert exterior ring
        exterior360 = [(x + 360.0 if x < 0 else x, y) for x, y in poly.exterior.coords]
        # convert any interior rings (holes)
        interiors360 = []
        for ring in poly.interiors:
            interiors360.append([(x + 360.0 if x < 0 else x, y) for x, y in ring.coords])
        return Polygon(exterior360, interiors360)    
    
    if isinstance(proxy_vals, xr.DataArray):
        lon_name = next((n for n in ["lon","longitude"] if n in proxy_vals.coords), None)
        lat_name = next((n for n in ["lat","latitude"] if n in proxy_vals.coords), None)

        if lon_name is None or lat_name is None:
            raise ValueError("Could not find longitude/latitude coordinates in DataArray.")

        lon_arr = proxy_vals[lon_name].values
        lat_arr = proxy_vals[lat_name].values

    if np.any(lon_arr < 0):
        lon_arr = np.where(lon_arr < 0, lon_arr + 360.0, lon_arr)

    # apply conversion to all polygons in the region
    polys360 = []
    for p in region.polygons:
        try:
            polys360.append(_poly_to_360(p))
        except Exception:
            # fallback: keep original polygon if conversion fails for any reason
            polys360.append(p)

    region_poly = unary_union(polys360)
    points = [Point(lon,lat) for lon, lat in zip(lon_arr,lat_arr)]
    mask = np.array([region_poly.contains(pt) for pt in points])

    if landsea_mask is not None:
        if not isinstance(landsea_mask, xr.DataArray):
            raise ValueError("landsea_mask must be an xarray.DataArray")
        
        lsm_lon_name = next((n for n in ["lon", "longitude"] if n in landsea_mask.coords), None)
        lsm_lat_name = next((n for n in ["lat", "latitude"] if n in landsea_mask.coords), None)

        if lsm_lon_name is None or lsm_lat_name is None:
            raise ValueError("Could not find longitude/latitude coordinates in landsea_mask DataArray.")

        lsm_lon = landsea_mask[lsm_lon_name].values
        lsm_lat = landsea_mask[lsm_lat_name].values 

        if np.any(lsm_lon < 0):
            lsm_lon = np.where(lsm_lon < 0, lsm_lon + 360.0, lsm_lon)

        lon_idx = np.argmin(np.abs(lsm_lon[:, None] - lon_arr[None, :]), axis=0)
        lat_idx = np.argmin(np.abs(lsm_lat[:, None] - lat_arr[None, :]), axis=0)

        land_vals = landsea_mask.values[lat_idx, lon_idx]

        land_mask = np.asarray(land_vals, dtype=bool)
        mask = mask & land_mask

    if isinstance(proxy_vals, xr.DataArray):
        mask_da = xr.DataArray(mask, coords=[proxy_vals.coords["points"]], dims=["points"])
        return proxy_vals.where(mask_da, drop=True)

    else:
        proxy_vals_out = proxy_vals[mask]
        proxy_lon_out = lon_arr[mask]
        proxy_lat_out = lat_arr[mask]

        return proxy_vals_out, proxy_lon_out, proxy_lat_out
    

