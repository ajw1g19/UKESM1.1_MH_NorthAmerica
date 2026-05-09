"""
Workaround for Cartopy + Shapely 2.x geometry errors.

This module patches shapely.geometry.MultiPolygon to handle two errors:

1. TypeError: 'GeometryCollection' object is not subscriptable
   - Caused by Shapely 2.x no longer supporting subscript on GeometryCollection
   - Solution: Extract Polygon objects from GeometryCollection before construction

2. ValueError: Sequences of multi-polygons are not valid arguments
   - Caused by Shapely 2.x rejecting MultiPolygon objects in input sequences
   - Note: MultiPolygon does NOT inherit from GeometryCollection in Shapely 2.x
   - Solution: Flatten MultiPolygon objects by extracting their constituent Polygons

Both errors occur during contourf operations with projection transformations
when Cartopy's _rings_to_multi_polygon() produces problematic geometry lists.

Issue: https://github.com/SciTools/cartopy/issues/2176
Remove this patch once Cartopy officially resolves the issue.
"""

import logging
from shapely.geometry import MultiPolygon, Polygon, GeometryCollection, LineString, Point

logger = logging.getLogger(__name__)

# Store original constructor
_original_multipolygon_new = MultiPolygon.__new__
_patch_applied = False


def _filtered_multipolygon_new(cls, polygons=None, **kwargs):
    """
    Patched MultiPolygon constructor that handles problematic geometry inputs.

    Handles two cases from Cartopy's projection transformations:

    1. GeometryCollection objects in the input list
       - Shapely 2.x no longer supports subscripting GeometryCollection
       - Solution: Extract Polygon objects from the GeometryCollection

    2. MultiPolygon objects in the input list
       - Shapely 2.x raises ValueError for sequences containing MultiPolygons
       - Note: MultiPolygon does NOT inherit from GeometryCollection in Shapely 2.x
       - Solution: Flatten by extracting constituent Polygon objects

    This patch extracts Polygon objects from problematic geometry types
    before passing to the original constructor.
    """
    if polygons is not None:
        # Check if any items need filtering:
        # - GeometryCollection (but not Polygon/MultiPolygon which are separate classes)
        # - MultiPolygon at top level (must be flattened)
        needs_filtering = any(
            (isinstance(item, GeometryCollection) and not isinstance(item, (Polygon, MultiPolygon)))
            or isinstance(item, MultiPolygon)
            for item in polygons
        )

        if needs_filtering:
            filtered_polygons = []

            for item in polygons:
                if isinstance(item, MultiPolygon):
                    # Flatten MultiPolygon to individual Polygons
                    filtered_polygons.extend(item.geoms)
                    logger.warning("MultiPolygon flattening patch applied")
                elif isinstance(item, GeometryCollection) and not isinstance(item, Polygon):
                    # Extract only Polygon/MultiPolygon from GeometryCollection
                    for geom in item.geoms:
                        if isinstance(geom, Polygon):
                            filtered_polygons.append(geom)
                        elif isinstance(geom, MultiPolygon):
                            filtered_polygons.extend(geom.geoms)
                        elif isinstance(geom, (LineString, Point)):
                            logger.debug(f"GeometryCollection filter: skipped {type(geom).__name__}")
                        else:
                            logger.debug(f"GeometryCollection filter: added {type(geom).__name__}")
                            filtered_polygons.append(geom)
                    # logger.warning("GeometryCollection patch applied")
                else:
                    # Keep Polygon and coordinate tuples as-is
                    filtered_polygons.append(item)

            polygons = filtered_polygons if filtered_polygons else polygons

    return _original_multipolygon_new(cls, polygons, **kwargs)


def apply_geometry_collection_fix():
    """
    Apply the geometry filter patch to shapely.MultiPolygon.

    This patches MultiPolygon.__new__ to handle:
    1. GeometryCollection objects (TypeError: not subscriptable)
    2. MultiPolygon objects in sequences (ValueError: not valid arguments)

    This should be called once at application startup, before any
    plotting operations that might trigger the error.

    Returns:
        bool: True if patch was applied, False if already applied
    """
    global _patch_applied

    if _patch_applied:
        logger.debug("GeometryCollection fix already applied")
        return False

    # Patch the __new__ method directly (not as staticmethod)
    MultiPolygon.__new__ = _filtered_multipolygon_new  # type: ignore[assignment]
    _patch_applied = True

    logger.debug(
        "Applied geometry filter patch for Cartopy/Shapely 2.x compatibility "
        "(handles GeometryCollection and MultiPolygon)"
    )
    return True


def remove_geometry_collection_fix():
    """
    Remove the patch and restore original MultiPolygon constructor.

    Useful for testing or when upgrading to a fixed version of Cartopy.
    """
    global _patch_applied

    if not _patch_applied:
        return False

    # Restore the original __new__ method
    MultiPolygon.__new__ = _original_multipolygon_new
    _patch_applied = False

    logger.info("Removed GeometryCollection filter patch")
    return True
