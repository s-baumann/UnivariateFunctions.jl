using Test
@testset "Date Tests" begin
    using UnivariateFunctions
    using Dates, TimeZones
    tol = 10*eps()

    # ========== Test unix2zdt ==========
    @testset "unix2zdt" begin
        # Test with integer input
        zdt = unix2zdt(0)
        @test zdt isa ZonedDateTime
        @test zdt == ZonedDateTime(1970, 1, 1, 0, 0, 0, tz"UTC")

        # Test with a known timestamp (2020-01-01 00:00:00 UTC = 1577836800)
        zdt_2020 = unix2zdt(1577836800)
        @test year(zdt_2020) == 2020
        @test month(zdt_2020) == 1
        @test day(zdt_2020) == 1

        # Test with different timezone
        zdt_ny = unix2zdt(0, tz"America/New_York")
        @test zdt_ny isa ZonedDateTime

        # Test with Real input (should floor)
        zdt_real = unix2zdt(1577836800.7)
        @test zdt_real == unix2zdt(1577836800)

        # Test with missing
        @test unix2zdt(missing) === missing

        # Test error for far future timestamps (likely microseconds mistake)
        @test_throws ErrorException unix2zdt(3405888001)
    end

    # ========== Test unix2dt ==========
    @testset "unix2dt" begin
        # Test with integer input
        dt = unix2dt(0)
        @test dt isa DateTime
        @test dt == DateTime(1970, 1, 1, 0, 0, 0)

        # Test with a known timestamp
        dt_2020 = unix2dt(1577836800)
        @test year(dt_2020) == 2020
        @test month(dt_2020) == 1
        @test day(dt_2020) == 1

        # Test with Real input (should floor)
        dt_real = unix2dt(1577836800.9)
        @test dt_real == unix2dt(1577836800)

        # Test with missing
        @test unix2dt(missing) === missing

        # Test error for far future timestamps
        @test_throws ErrorException unix2dt(3405888001)
    end

    # ========== Test unix2d ==========
    @testset "unix2d" begin
        # Test with integer input
        d = unix2d(0)
        @test d isa Date
        @test d == Date(1970, 1, 1)

        # Test with a known timestamp
        d_2020 = unix2d(1577836800)
        @test d_2020 == Date(2020, 1, 1)

        # Test with Real input (should floor)
        d_real = unix2d(1577836800.9)
        @test d_real == unix2d(1577836800)

        # Test with missing
        @test unix2d(missing) === missing

        # Test error for far future timestamps
        @test_throws ErrorException unix2d(3405888001)
    end

    # ========== Test zdt2unix ==========
    @testset "zdt2unix" begin
        # Test with ZonedDateTime
        zdt = ZonedDateTime(1970, 1, 1, 0, 0, 0, tz"UTC")
        @test zdt2unix(zdt) == 0

        zdt_2020 = ZonedDateTime(2020, 1, 1, 0, 0, 0, tz"UTC")
        @test zdt2unix(zdt_2020) == 1577836800

        # Test with DateTime
        dt = DateTime(1970, 1, 1, 0, 0, 0)
        @test zdt2unix(dt) == 0

        dt_2020 = DateTime(2020, 1, 1, 0, 0, 0)
        @test zdt2unix(dt_2020) == 1577836800

        # Test with Date
        d = Date(1970, 1, 1)
        @test zdt2unix(d) == 0

        d_2020 = Date(2020, 1, 1)
        @test zdt2unix(d_2020) == 1577836800

        # Test with missing
        @test zdt2unix(missing) === missing

        # Test roundtrip
        original_unix = 1600000000
        @test zdt2unix(unix2zdt(original_unix)) == original_unix
        @test zdt2unix(unix2dt(original_unix)) == original_unix
        @test zdt2unix(unix2d(original_unix)) == 1599955200  # Date loses time precision
    end

    # ========== Test seconds_between ==========
    @testset "seconds_between" begin
        # Test with Dates
        d1 = Date(2020, 1, 1)
        d2 = Date(2020, 1, 2)
        @test seconds_between(d2, d1) == 86400  # 1 day = 86400 seconds

        # Test with DateTimes
        dt1 = DateTime(2020, 1, 1, 0, 0, 0)
        dt2 = DateTime(2020, 1, 1, 1, 0, 0)
        @test seconds_between(dt2, dt1) == 3600  # 1 hour = 3600 seconds

        # Test with ZonedDateTimes
        zdt1 = ZonedDateTime(2020, 1, 1, 0, 0, 0, tz"UTC")
        zdt2 = ZonedDateTime(2020, 1, 1, 0, 1, 0, tz"UTC")
        @test seconds_between(zdt2, zdt1) == 60  # 1 minute = 60 seconds

        # Test negative (b > a)
        @test seconds_between(d1, d2) == -86400

        # Test mixed types
        @test seconds_between(dt2, d1) isa Real
        @test seconds_between(zdt2, dt1) isa Real
    end

    # ========== Test days_between ==========
    @testset "days_between" begin
        d1 = Date(2020, 1, 1)
        d2 = Date(2020, 1, 11)
        @test days_between(d2, d1) == 10.0

        # Test with partial days using DateTime
        dt1 = DateTime(2020, 1, 1, 0, 0, 0)
        dt2 = DateTime(2020, 1, 1, 12, 0, 0)
        @test days_between(dt2, dt1) == 0.5

        # Test negative
        @test days_between(d1, d2) == -10.0
    end

    # ========== Test years_between ==========
    @testset "years_between" begin
        d1 = Date(2020, 1, 1)
        d2 = Date(2021, 1, 1)
        # 366 days in 2020 (leap year), so slightly more than 1 year
        yrs = years_between(d2, d1)
        @test yrs > 1.0
        @test yrs < 1.01  # Should be close to 1

        # Test longer period
        d3 = Date(2030, 1, 1)
        yrs_10 = years_between(d3, d1)
        @test yrs_10 > 9.9
        @test yrs_10 < 10.1

        # Test negative
        @test years_between(d1, d2) < 0
    end

    # ========== Test period_length ==========
    @testset "period_length" begin
        # Test with Year
        @test period_length(Year(1)) ≈ 1.0 atol=0.01

        # Test with Month (approximately 1/12 of a year)
        month_len = period_length(Month(1))
        @test month_len > 0.08
        @test month_len < 0.09

        # Test with Day
        day_len = period_length(Day(1))
        @test day_len ≈ 1/365.25 atol=0.0001

        # Test with Week
        week_len = period_length(Week(1))
        @test week_len ≈ 7/365.25 atol=0.001

        # Test multiple periods
        @test period_length(Year(5)) ≈ 5.0 atol=0.05
        @test period_length(Month(12)) ≈ 1.0 atol=0.01
    end

    # ========== Test years_from_global_base_date ==========
    @testset "years_from_global_base_date" begin
        # Test at the global base date (1970-01-01)
        @test years_from_global_base_date(Date(1970, 1, 1)) == 0.0

        # Test one year later
        yrs_1971 = years_from_global_base_date(Date(1971, 1, 1))
        @test yrs_1971 > 0.99
        @test yrs_1971 < 1.01

        # Test with DateTime
        yrs_dt = years_from_global_base_date(DateTime(1970, 1, 1, 12, 0, 0))
        @test yrs_dt > 0.0
        @test yrs_dt < 0.01  # Half a day

        # Test with ZonedDateTime
        yrs_zdt = years_from_global_base_date(ZonedDateTime(1970, 1, 1, 0, 0, 0, tz"UTC"))
        @test yrs_zdt == 0.0

        # Test before base date (negative)
        yrs_before = years_from_global_base_date(Date(1969, 1, 1))
        @test yrs_before < 0.0
        @test yrs_before > -1.1

        # Test far in future
        yrs_2050 = years_from_global_base_date(Date(2050, 1, 1))
        @test yrs_2050 > 79.9
        @test yrs_2050 < 80.1
    end

    # ========== Test PE_Function with dates (original tests) ==========
    @testset "PE_Function with dates" begin
        today = Date(2000,1,1)
        pe_func = PE_Function(1.0,2.0,today, 3)
        @test abs(pe_func.base_ - years_from_global_base_date(today)) < tol

        date_in_2020 = Date(2000,2,3)
        pe_func2 = PE_Function(1.0,0.00000000002,date_in_2020, 2)
        @test abs(pe_func2.base_ - years_from_global_base_date(date_in_2020)) < tol
        @test abs(evaluate(pe_func, date_in_2020) - evaluate(pe_func, years_from_global_base_date(date_in_2020))) < tol

        # Sum of functions
        sum_func = Sum_Of_Functions([pe_func, PE_Function(2.0,0.00000000000025,today, 1)])
        @test abs(evaluate(sum_func, date_in_2020) - evaluate(sum_func, years_from_global_base_date(date_in_2020))) < tol

        # Left and right integrals
        l_int = left_integral(pe_func, today)
        @test (evaluate(l_int, date_in_2020) - evaluate_integral(pe_func, today, date_in_2020)) < tol
        r_int = right_integral(pe_func, date_in_2020)
        @test (evaluate(r_int, today) - evaluate_integral(pe_func, today, date_in_2020)) < tol

        # With DateTimes
        today_time = DateTime(2000,1,1, 10, 10, 1)
        pe_func = PE_Function(1.0,0.00000000000020,today_time, 1)

        # Left and right integrals
        later_today = DateTime(2000,1,1, 10, 14, 1)
        l_int = left_integral(pe_func, today)
        @test (evaluate(l_int, later_today) - evaluate_integral(pe_func, today, later_today)) < tol
        r_int = right_integral(pe_func, later_today)
        @test (evaluate(r_int, today) - evaluate_integral(pe_func, today, later_today)) < tol
    end

    # ========== Test edge cases ==========
    @testset "Edge cases" begin
        # Test with leap year boundaries
        leap_day = Date(2020, 2, 29)
        @test years_from_global_base_date(leap_day) isa Real

        # Test consistency across types
        d = Date(2020, 6, 15)
        dt = DateTime(2020, 6, 15, 0, 0, 0)
        zdt = ZonedDateTime(2020, 6, 15, 0, 0, 0, tz"UTC")

        @test zdt2unix(d) == zdt2unix(dt)
        @test zdt2unix(dt) == zdt2unix(zdt)

        # Test that unix2* and zdt2unix are inverses (within precision)
        unix_time = 1592179200  # 2020-06-15 00:00:00 UTC
        @test zdt2unix(unix2zdt(unix_time)) == unix_time
        @test zdt2unix(unix2dt(unix_time)) == unix_time

        # Test seconds_between with same date
        @test seconds_between(d, d) == 0
        @test days_between(d, d) == 0.0
        @test years_between(d, d) == 0.0
    end

    println("Date conversion tests passed.")
end
