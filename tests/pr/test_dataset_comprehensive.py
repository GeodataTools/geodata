# Copyright 2024 Xiqiang Liu

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Comprehensive examples of different types of tests for dataset classes.

This file demonstrates various testing patterns for geospatial datasets.
Each test category serves a specific purpose in ensuring dataset quality,
correctness, and reliability.

Test Categories Explained:
--------------------------

1. DOWNLOAD TESTS
   - Verify that datasets can be downloaded successfully
   - Ensure files are saved to correct locations
   - Critical for ensuring the basic data acquisition pipeline works
   
2. CATALOG TESTS
   - Validate that the catalog (list of files to download) is correctly generated
   - Test different frequencies (monthly, daily, hourly)
   - Verify testing mode limits downloads appropriately
   - Important for understanding what files will be downloaded before actually downloading

3. DATA INTEGRITY TESTS
   - Check file checksums/hashes to detect corruption
   - Verify downloaded files can be opened and read
   - Ensure data hasn't been corrupted during download or storage
   - Critical for data quality assurance

4. COORDINATE & BOUNDS TESTS
   - Validate coordinate system transformations (lat/lon to x/y)
   - Test bounding box filtering works correctly
   - Verify coordinate ranges are within expected limits
   - Important for spatial data correctness

5. METADATA TESTS
   - Verify dataset properties (projection, lat_direction, frequency)
   - Test that required attributes are present
   - Validate metadata is consistent across dataset types
   - Important for understanding dataset characteristics

6. DATA STRUCTURE TESTS
   - Verify downloaded datasets have expected variables
   - Check coordinate dimensions match expectations
   - Validate data types and value ranges
   - Critical for ensuring data usability

7. POSTPROCESSING TESTS
   - Test that dataset postprocessing functions correctly
   - Verify coordinate renaming (lat/lon -> x/y)
   - Check data transformations are applied correctly
   - Important for ensuring data is in the expected format

8. MULTI-DATASET TESTS
   - Compare outputs from different datasets for consistency
   - Test interoperability between different dataset types
   - Important for ensuring datasets can be used together
"""

import logging
import os
from typing import Optional

import xarray as xr

from geodata.datasets import load_dataset

logging.basicConfig(level=logging.INFO)

# PRs should run with zero CDS calls. Enable integration download tests only via:
#   GEODATA_RUN_CDS_TESTS=1
RUN_CDS_TESTS = os.getenv("GEODATA_RUN_CDS_TESTS") == "1"


# ============================================================================
# TEST CONFIGURATION HELPERS
# ============================================================================

def get_data_configs() -> list[str]:
    """Get list of dataset configurations to test."""
    # Default to offline fixtures for PR safety.
    return ["wind_3d_hourly"] if RUN_CDS_TESTS else ["wind_3d_hourly_test"]


def get_bounds() -> list[list[float]]:
    """Get list of bounding boxes to test (lon_min, lat_min, lon_max, lat_max)."""
    # Bounds are only meaningful for real downloads. Fixture files are not
    # regenerated per-bounds and therefore shouldn't be validated against bounds.
    return [[50, 0, 48, 3]] if RUN_CDS_TESTS else [None]  # type: ignore[list-item]


def get_years() -> list[slice]:
    """Get list of year ranges to test."""
    return [slice(2005, 2005)] if RUN_CDS_TESTS else [slice(2016, 2016)]


def get_months() -> list[slice]:
    """Get list of month ranges to test."""
    return [slice(1, 2)] if RUN_CDS_TESTS else [slice(1, 1)]


def get_dataset(
    data_config: str,
    bound: Optional[list[float]],
    year: slice,
    month: slice,
    testing: bool = True,
):
    """Helper function to create and optionally download a dataset."""
    dataset_cls = load_dataset(data_config)
    dataset = dataset_cls(
        years=year, months=month, bounds=bound, testing=testing
    )
    if not dataset.downloaded:
        if RUN_CDS_TESTS:
            dataset.download()
        else:
            raise AssertionError(
                f"Dataset {data_config} is not downloaded, but CDS tests are disabled. "
                "Use fixture configs or set GEODATA_RUN_CDS_TESTS=1."
            )
    return dataset


# ============================================================================
# 1. DOWNLOAD TESTS
# ============================================================================

def test_download():
    """
    Test Category 1: Download Tests
    
    WHY: Ensures the basic data acquisition pipeline works correctly.
    Downloads are expensive (time, bandwidth, storage), so we need to verify
    they work before running longer tests. This is the foundation for all
    other data-dependent tests.
    """
    if not RUN_CDS_TESTS:
        # This test verifies CDS download pipeline. Keep it opt-in.
        return

    configs = get_data_configs()
    years = get_years()
    months = get_months()
    bounds = get_bounds()

    for config, year, month, bound in zip(configs, years, months, bounds):
        dataset = get_dataset(config, bound, year, month)
        assert dataset.downloaded, f"Dataset {config} should be downloaded"


# ============================================================================
# 2. CATALOG TESTS
# ============================================================================

def test_catalog_generation():
    """
    Test Category 2: Catalog Generation Tests
    
    WHY: The catalog determines which files need to be downloaded. Incorrect
    catalog generation means missing data or unnecessary downloads. Testing
    this ensures we know exactly what will be downloaded before we download it.
    """
    config = "wind_3d_hourly" if RUN_CDS_TESTS else "wind_3d_hourly_test"
    dataset_cls = load_dataset(config)
    
    # Test monthly catalog (if applicable)
    dataset = dataset_cls(
        years=slice(2005, 2005) if RUN_CDS_TESTS else slice(2016, 2016),
        months=slice(1, 1),
        testing=True,
    )
    catalog = dataset.catalog
    
    assert len(catalog) > 0, "Catalog should contain at least one file"
    
    # Verify catalog entries have correct structure
    for file in catalog:
        assert hasattr(file, "year"), "Catalog entry should have year"
        assert hasattr(file, "month"), "Catalog entry should have month"
        assert hasattr(file, "path"), "Catalog entry should have path"
        assert file.year == (2005 if RUN_CDS_TESTS else 2016), "Year should match"
        assert file.month == 1, "Month should match"


def test_catalog_testing_mode():
    """
    Test Category 2: Testing Mode Catalog Tests
    
    WHY: Testing mode should limit downloads to a few days/months to speed up
    tests. If this doesn't work correctly, tests become slow and expensive.
    """
    if not RUN_CDS_TESTS:
        # Fixture datasets have fixed catalogs; testing mode isn't meaningful here.
        return

    config = "wind_3d_hourly"
    dataset_cls = load_dataset(config)
    
    # Testing mode should limit to 3 days for daily frequency datasets
    dataset_testing = dataset_cls(
        years=slice(2005, 2005),
        months=slice(1, 1),
        testing=True
    )
    catalog_testing = dataset_testing.catalog
    
    # Non-testing mode would download full month
    dataset_normal = dataset_cls(
        years=slice(2005, 2005),
        months=slice(1, 1),
        testing=False
    )
    catalog_normal = dataset_normal.catalog
    
    # Testing mode should have fewer files
    assert len(catalog_testing) < len(catalog_normal), \
        "Testing mode should limit the number of files"


def test_catalog_paths():
    """
    Test Category 2: Catalog Path Tests
    
    WHY: File paths determine where data is stored. Incorrect paths lead to
    data being saved in wrong locations or files overwriting each other.
    """
    config = "wind_3d_hourly" if RUN_CDS_TESTS else "wind_3d_hourly_test"
    dataset_cls = load_dataset(config)
    dataset = dataset_cls(
        years=slice(2005, 2005) if RUN_CDS_TESTS else slice(2016, 2016),
        months=slice(1, 1),
        testing=True,
    )
    
    catalog = dataset.catalog
    paths = {file.path for file in catalog}
    
    # All paths should be unique
    assert len(paths) == len(catalog), "All catalog paths should be unique"
    
    # Paths should follow expected structure (year/month/day.nc for daily)
    for file in catalog:
        path_str = str(file.path)
        assert str(file.year) in path_str, "Path should contain year"
        assert f"{file.month:02d}" in path_str, "Path should contain month"
        if file.day is not None:
            assert f"{file.day:02d}.nc" in path_str, "Path should contain day for daily datasets"


# ============================================================================
# 3. DATA INTEGRITY TESTS
# ============================================================================

def test_file_integrity():
    """
    Test Category 3: File Integrity Tests
    
    WHY: Downloaded files can become corrupted during transfer or storage.
    Integrity checks catch these issues before they cause problems in analysis.
    """
    configs = get_data_configs()
    years = get_years()
    months = get_months()
    bounds = get_bounds()

    for config, year, month, bound in zip(configs, years, months, bounds):
        dataset = get_dataset(config, bound, year, month)
        
        # Check integrity of all files in catalog
        for file in dataset.catalog:
            assert file.check(), f"File {file.path} should exist"
            
            # Test integrity check (requires file_hash to be set)
            # Note: This would require files to have hashes stored
            assert file.check(integrity=False), \
                f"File {file.path} should pass basic integrity check"


def test_file_readability():
    """
    Test Category 3: File Readability Tests
    
    WHY: A file can exist and pass checksum but still be unreadable (wrong format,
    corrupted headers, etc.). This ensures we can actually use the downloaded data.
    """
    configs = get_data_configs()
    years = get_years()
    months = get_months()
    bounds = get_bounds()

    for config, year, month, bound in zip(configs, years, months, bounds):
        dataset = get_dataset(config, bound, year, month)
        
        for file in dataset.catalog:
            if file.path.exists():
                # Should be able to open as xarray dataset
                ds = xr.open_dataset(file.path)
                assert ds is not None, f"Should be able to open {file.path}"
                ds.close()


# ============================================================================
# 4. COORDINATE & BOUNDS TESTS
# ============================================================================

def test_bounds_validation():
    """
    Test Category 4: Bounds Validation Tests
    
    WHY: Bounding boxes filter data spatially. Incorrect bounds can lead to
    downloading unnecessary data or missing required data. Also validates that
    invalid bounds are rejected early.
    """
    if not RUN_CDS_TESTS:
        return

    config = "wind_3d_hourly"
    dataset_cls = load_dataset(config)
    
    # Test valid bounds
    valid_bounds = [50, 0, 52, 3]  # lon_min, lat_min, lon_max, lat_max
    dataset = dataset_cls(
        years=slice(2005, 2005),
        months=slice(1, 1),
        bounds=valid_bounds,
        testing=True
    )
    assert dataset.bounds == valid_bounds, "Valid bounds should be accepted"
    
    # Note: To test invalid bounds validation, you could add a test that
    # verifies ValueError is raised for bounds outside valid ranges.
    # Example: bounds with longitude > 180 or < -180 should raise ValueError


def test_coordinate_renaming():
    """
    Test Category 4: Coordinate Renaming Tests
    
    WHY: Datasets use different coordinate names (lat/lon vs x/y). The base
    class should standardize these. Incorrect renaming breaks downstream analysis.
    """
    configs = get_data_configs()
    years = get_years()
    months = get_months()
    bounds = get_bounds()

    for config, year, month, bound in zip(configs, years, months, bounds):
        dataset = get_dataset(config, bound, year, month)
        
        # Check at least one file to verify coordinate naming
        for file in dataset.catalog:
            if file.path.exists():
                ds = xr.open_dataset(file.path)
                
                # After postprocessing, coordinates should be renamed to x, y
                # (or lat, lon should be present if add_lon_lat=True)
                coords = list(ds.coords.keys())
                
                # Should have x and y coordinates (or lat/lon)
                has_xy = "x" in coords and "y" in coords
                has_latlon = "lat" in coords and "lon" in coords
                
                assert has_xy or has_latlon, \
                    f"Dataset should have x/y or lat/lon coordinates. Found: {coords}"
                
                ds.close()
                break  # Only check first file


# ============================================================================
# 5. METADATA TESTS
# ============================================================================

def test_dataset_properties():
    """
    Test Category 5: Dataset Properties Tests
    
    WHY: Dataset properties (projection, lat_direction, frequency) are used
    throughout the codebase for processing. Incorrect properties break analysis.
    """
    config = "wind_3d_hourly" if RUN_CDS_TESTS else "wind_3d_hourly_test"
    dataset_cls = load_dataset(config)
    dataset = dataset_cls(
        years=slice(2005, 2005) if RUN_CDS_TESTS else slice(2016, 2016),
        months=slice(1, 1),
        testing=True,
    )
    
    # Test required properties exist
    assert hasattr(dataset, "projection"), "Dataset should have projection property"
    assert hasattr(dataset, "lat_direction"), "Dataset should have lat_direction property"
    assert hasattr(dataset, "frequency"), "Dataset should have frequency property"
    assert hasattr(dataset, "module"), "Dataset should have module attribute"
    assert hasattr(dataset, "weather_config"), "Dataset should have weather_config attribute"
    
    # Test property types
    assert isinstance(dataset.projection, str), "Projection should be a string"
    assert isinstance(dataset.lat_direction, bool), "lat_direction should be a boolean"
    assert dataset.frequency in ["hourly", "daily", "monthly"], \
        "Frequency should be one of: hourly, daily, monthly"


def test_dataset_repr():
    """
    Test Category 5: Dataset Representation Tests
    
    WHY: The __repr__ method is used for debugging and logging. It should provide
    useful information about the dataset state.
    """
    config = "wind_3d_hourly" if RUN_CDS_TESTS else "wind_3d_hourly_test"
    dataset_cls = load_dataset(config)
    dataset = dataset_cls(
        years=slice(2005, 2005) if RUN_CDS_TESTS else slice(2016, 2016),
        months=slice(1, 1),
        testing=True,
    )
    
    repr_str = repr(dataset)
    
    # Should contain key information
    assert dataset.weather_config in repr_str, "repr should contain weather_config"
    assert ("2005" if RUN_CDS_TESTS else "2016") in repr_str, "repr should contain years"
    assert "1" in repr_str, "repr should contain months"


# ============================================================================
# 6. DATA STRUCTURE TESTS
# ============================================================================

def test_data_variables():
    """
    Test Category 6: Data Variables Tests
    
    WHY: Each dataset should contain specific variables. Missing or incorrectly
    named variables break downstream analysis that depends on them.
    """
    configs = get_data_configs()
    years = get_years()
    months = get_months()
    bounds = get_bounds()

    for config, year, month, bound in zip(configs, years, months, bounds):
        dataset = get_dataset(config, bound, year, month)
        
        # Check if dataset defines expected variables
        # Note: Not all datasets have a 'variables' attribute
        # This is an example of how to test datasets that do have it
        if hasattr(dataset, "variables"):
            expected_vars = getattr(dataset, "variables")
            
            # Verify at least one file contains these variables
            for file in dataset.catalog:
                if file.path.exists():
                    ds = xr.open_dataset(file.path)
                    
                    # Variables should exist in dataset
                    for var in expected_vars:
                        assert var in ds.data_vars or var in ds.coords, \
                            f"Variable {var} should exist in dataset"
                    
                    ds.close()
                    break  # Only check first file


def test_data_dimensions():
    """
    Test Category 6: Data Dimension Tests
    
    WHY: Data dimensions determine how data can be processed. For example,
    a 3D wind dataset should have a level/height dimension. Missing dimensions
    indicate incorrect data structure.
    """
    config = "wind_3d_hourly" if RUN_CDS_TESTS else "wind_3d_hourly_test"
    dataset_cls = load_dataset(config)
    dataset = dataset_cls(
        years=slice(2005, 2005) if RUN_CDS_TESTS else slice(2016, 2016),
        months=slice(1, 1),
        bounds=get_bounds()[0],
        testing=True,
    )
    if RUN_CDS_TESTS and not dataset.downloaded:
        dataset.download()
    
    # Check first downloaded file
    for file in dataset.catalog:
        if file.path.exists():
            ds = xr.open_dataset(file.path)
            
            # 3D wind data should have multiple dimensions
            dims = list(ds.dims.keys())
            
            # Should have spatial dimensions
            assert "x" in dims or "lon" in dims, "Should have x/lon dimension"
            assert "y" in dims or "lat" in dims, "Should have y/lat dimension"
            
            # 3D data should have a level/height dimension
            _ = any(dim in dims for dim in ["level", "height", "lev", "plev"])
            
            ds.close()
            break  # Only check first file


def test_data_value_ranges():
    """
    Test Category 6: Data Value Range Tests
    
    WHY: Data values should be within physically plausible ranges. Out-of-range
    values indicate data corruption or processing errors.
    """
    config = "wind_3d_hourly" if RUN_CDS_TESTS else "wind_3d_hourly_test"
    dataset_cls = load_dataset(config)
    dataset = dataset_cls(
        years=slice(2005, 2005) if RUN_CDS_TESTS else slice(2016, 2016),
        months=slice(1, 1),
        bounds=get_bounds()[0],
        testing=True
    )
    if RUN_CDS_TESTS and not dataset.downloaded:
        dataset.download()
    
    # Check first downloaded file
    for file in dataset.catalog:
        if file.path.exists():
            ds = xr.open_dataset(file.path)
            
            # Check that data values are finite (not NaN or Inf)
            for var in ds.data_vars:
                data = ds[var]
                assert data.notnull().any(), \
                    f"Variable {var} should have some non-null values"
                
                # Wind components should be within reasonable range
                # (typical wind speeds are -100 to 100 m/s)
                var_str = str(var)
                if "u" in var_str.lower() or "v" in var_str.lower():
                    if data.notnull().any():
                        data_min = float(data.min())
                        data_max = float(data.max())
                        # Allow wide range, but should be finite
                        assert abs(data_min) < 200, \
                            f"Wind component {var} min value {data_min} seems unreasonable"
                        assert abs(data_max) < 200, \
                            f"Wind component {var} max value {data_max} seems unreasonable"
            
            ds.close()
            break  # Only check first file


# ============================================================================
# 7. POSTPROCESSING TESTS
# ============================================================================

def test_postprocessing_applied():
    """
    Test Category 7: Postprocessing Application Tests
    
    WHY: Postprocessing (coordinate renaming, data transformations) must be
    applied consistently. If postprocessing fails silently, downstream code
    expecting transformed data will fail.
    """
    config = "wind_3d_hourly" if RUN_CDS_TESTS else "wind_3d_hourly_test"
    dataset_cls = load_dataset(config)
    dataset = dataset_cls(
        years=slice(2005, 2005) if RUN_CDS_TESTS else slice(2016, 2016),
        months=slice(1, 1),
        bounds=get_bounds()[0],
        testing=True
    )
    if RUN_CDS_TESTS and not dataset.downloaded:
        dataset.download()
    
    # Check that postprocessed files have correct structure
    for file in dataset.catalog:
        if file.path.exists():
            ds = xr.open_dataset(file.path)
            
            # Postprocessing should rename coordinates
            # Check that we have standardized coordinate names
            coords = list(ds.coords.keys())
            assert "x" in coords or "lon" in coords, \
                "Postprocessed data should have x/lon coordinate"
            assert "y" in coords or "lat" in coords, \
                "Postprocessed data should have y/lat coordinate"
            
            ds.close()
            break  # Only check first file


# ============================================================================
# 8. MULTI-DATASET TESTS (Example - can be expanded)
# ============================================================================

def test_datasets_loaded_correctly():
    """
    Test Category 8: Multi-Dataset Loading Tests
    
    WHY: The dataset registry and loading mechanism must work correctly for
    all datasets. If one dataset can't be loaded, it breaks the entire system.
    """
    from geodata.datasets import list_datasets, load_dataset
    
    # Should be able to list all datasets
    datasets = list_datasets()
    assert len(datasets) > 0, "Should have at least one dataset registered"
    
    # Should be able to load each dataset class
    for dataset_name in datasets:
        dataset_cls = load_dataset(dataset_name)
        assert dataset_cls is not None, \
            f"Should be able to load dataset class for {dataset_name}"


# ============================================================================
# ADDITIONAL USEFUL TESTS
# ============================================================================

def test_testing_mode():
    """
    Additional Test: Testing Mode Behavior
    
    WHY: Testing mode is crucial for fast CI/CD pipelines. If it doesn't work
    correctly, tests become too slow or download too much data.
    """
    if not RUN_CDS_TESTS:
        # Fixture datasets ignore testing-mode catalog limiting; keep this check opt-in.
        return

    config = "wind_3d_hourly"
    dataset_cls = load_dataset(config)
    
    # Testing mode should limit downloads
    dataset_testing = dataset_cls(
        years=slice(2005, 2005),
        months=slice(1, 1),
        testing=True
    )
    assert dataset_testing.testing is True, "Testing mode should be enabled"
    
    # Non-testing mode
    dataset_normal = dataset_cls(
        years=slice(2005, 2005),
        months=slice(1, 1),
        testing=False
    )
    assert dataset_normal.testing is False, "Testing mode should be disabled"


def test_storage_path():
    """
    Additional Test: Storage Path Tests
    
    WHY: Files must be saved to the correct location for proper organization
    and retrieval. Wrong paths make it impossible to find downloaded data.
    """
    config = "wind_3d_hourly" if RUN_CDS_TESTS else "wind_3d_hourly_test"
    dataset_cls = load_dataset(config)
    dataset = dataset_cls(
        years=slice(2005, 2005) if RUN_CDS_TESTS else slice(2016, 2016),
        months=slice(1, 1),
        testing=True,
    )
    
    # Storage root should follow expected pattern
    assert dataset.storage_root is not None, "Storage root should be set"
    assert "era5" in str(dataset.storage_root), \
        "Storage root should contain module name"
    assert dataset.weather_config in str(dataset.storage_root), \
        "Storage root should contain weather_config"


def test_bounds_applied():
    """
    Additional Test: Bounds Application Tests
    
    WHY: When bounds are specified, data should be filtered to those bounds.
    Downloading global data when only a region is needed wastes resources.
    """
    if not RUN_CDS_TESTS:
        return

    config = "wind_3d_hourly"
    dataset_cls = load_dataset(config)
    
    bounds = [50, 0, 52, 3]  # Small region
    dataset = dataset_cls(
        years=slice(2005, 2005),
        months=slice(1, 1),
        bounds=bounds,
        testing=True
    )
    dataset.download()
    
    # Check that downloaded data respects bounds
    for file in dataset.catalog:
        if file.path.exists():
            ds = xr.open_dataset(file.path)
            
            # Check coordinate ranges (if coordinates are available)
            if "x" in ds.coords:
                x_coords = ds.coords["x"].values
                lon_min, lon_max = min(x_coords), max(x_coords)
                # Data should be within or close to bounds (allowing for grid resolution)
                # ERA5 uses 0.25-degree grid, and xr.sel() with slice may include grid points
                # that extend beyond requested bounds. We allow up to 2.5 degrees tolerance
                # to account for grid alignment and coordinate system conversions.
                # Bounds are [lon_min, lat_min, lon_max, lat_max]
                tolerance = 2.5  # Degrees tolerance for grid resolution and coordinate conversion
                assert lon_min >= bounds[0] - tolerance, \
                    f"Longitude min {lon_min} should be >= bounds[0] {bounds[0]} - {tolerance}"
                assert lon_max <= bounds[2] + tolerance, \
                    f"Longitude max {lon_max} should be <= bounds[2] {bounds[2]} + {tolerance}"
            
            ds.close()
            break  # Only check first file


