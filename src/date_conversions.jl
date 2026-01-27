"""
    unix2zdt(seconds_since_1970, tzz=tz"UTC")

Convert Unix timestamp (seconds since 1970-01-01 00:00:00 UTC) to a `ZonedDateTime`.

### Inputs
* `seconds_since_1970` - Unix timestamp as Integer, Real, or missing
* `tzz` - Timezone (default: UTC)

### Returns
* A `ZonedDateTime`, or `missing` if input is missing.

### Example
```julia
unix2zdt(1577836800)  # 2020-01-01T00:00:00+00:00
unix2zdt(0, tz"America/New_York")  # 1969-12-31T19:00:00-05:00
```
"""
function unix2zdt(seconds_since_1970::Integer, tzz = tz"UTC")
    if seconds_since_1970 > 3405888000
        error("This seems wrong unless you are far into the future. Are you sure you are not inputting microseconds since 1970?")
    end
    return ZonedDateTime(unix2datetime(seconds_since_1970), tzz)
end

"""
    unix2dt(seconds_since_1970)

Convert Unix timestamp (seconds since 1970-01-01 00:00:00 UTC) to a `DateTime`.

### Inputs
* `seconds_since_1970` - Unix timestamp as Integer, Real, or missing

### Returns
* A `DateTime`, or `missing` if input is missing.

### Example
```julia
unix2dt(1577836800)  # 2020-01-01T00:00:00
```
"""
function unix2dt(seconds_since_1970::Integer)
    if seconds_since_1970 > 3405888000
        error("This seems wrong unless you are far into the future. Are you sure you are not inputting microseconds since 1970?")
    end
    return unix2datetime(seconds_since_1970)
end

"""
    unix2d(seconds_since_1970)

Convert Unix timestamp (seconds since 1970-01-01 00:00:00 UTC) to a `Date`.

### Inputs
* `seconds_since_1970` - Unix timestamp as Integer, Real, or missing

### Returns
* A `Date`, or `missing` if input is missing.

### Example
```julia
unix2d(1577836800)  # 2020-01-01
```
"""
function unix2d(seconds_since_1970::Integer)
    if seconds_since_1970 > 3405888000
        error("This seems wrong unless you are far into the future. Are you sure you are not inputting microseconds since 1970?")
    end
    return Date(unix2datetime(seconds_since_1970))
end
function unix2d(seconds_since_1970::Missing)
    return missing
end
function unix2dt(seconds_since_1970::Missing)
    return missing
end
function unix2zdt(seconds_since_1970::Missing)
    return missing
end
function unix2zdt(seconds_since_1970::Real, tzz = tz"UTC")
    return unix2zdt(Int(floor(seconds_since_1970)), tzz)
end
function unix2dt(seconds_since_1970::Real)
    return unix2dt(Int(floor(seconds_since_1970)))
end
function unix2d(seconds_since_1970::Real)
    return unix2d(Int(floor(seconds_since_1970)))
end

"""
    zdt2unix(d)

Convert a `ZonedDateTime`, `DateTime`, or `Date` to Unix timestamp (seconds since 1970-01-01 00:00:00 UTC).

### Inputs
* `d` - A `ZonedDateTime`, `DateTime`, `Date`, or `missing`

### Returns
* An `Int` representing seconds since Unix epoch, or `missing` if input is missing.

### Example
```julia
zdt2unix(Date(2020, 1, 1))  # 1577836800
zdt2unix(DateTime(2020, 1, 1, 12, 0, 0))  # 1577880000
```
"""
function zdt2unix(zdt::ZonedDateTime)
    return Int(floor(datetime2unix(zdt.utc_datetime)))
end
function zdt2unix(dt::DateTime)
    return Int(floor(datetime2unix(dt)))
end
function zdt2unix(d::Date)
    return Int(floor(datetime2unix(DateTime(d))))
end
function zdt2unix(d::Missing)
    return missing
end

#########

"""
    seconds_between(a::Union{ZonedDateTime,DateTime,Date}, b::Union{ZonedDateTime,DateTime,Date})
The number of seconds between two dates. This is the delta between two unixtimes.

### Inputs
* `a` - The end date
* `b` - The start date.
### Returns
* A `Real`.
"""
function seconds_between(a::Union{ZonedDateTime,DateTime,Date}, b::Union{ZonedDateTime,DateTime,Date})
    return zdt2unix(a) - zdt2unix(b)
end

"""
    days_between(a, b)

The number of days between two dates. Computed as `seconds_between(a, b) / 86400`.

### Inputs
* `a` - The end date (ZonedDateTime, DateTime, or Date)
* `b` - The start date (ZonedDateTime, DateTime, or Date)

### Returns
* A `Real` representing the number of days (can be fractional).

### Example
```julia
days_between(Date(2020, 1, 11), Date(2020, 1, 1))  # 10.0
days_between(DateTime(2020, 1, 1, 12, 0, 0), DateTime(2020, 1, 1, 0, 0, 0))  # 0.5
```
"""
function days_between(a::Union{ZonedDateTime,DateTime,Date}, b::Union{ZonedDateTime,DateTime,Date})
    return (seconds_between(a, b) / 86400)
end

"""
    years_between(a, b)

The number of years between two dates. Computed as `seconds_between(a, b) / (365.25 * 86400)`.

### Inputs
* `a` - The end date (ZonedDateTime, DateTime, or Date)
* `b` - The start date (ZonedDateTime, DateTime, or Date)

### Returns
* A `Real` representing the number of years (can be fractional).

### Example
```julia
years_between(Date(2021, 1, 1), Date(2020, 1, 1))  # approximately 1.0
years_between(Date(2020, 7, 1), Date(2020, 1, 1))  # approximately 0.5
```
"""
function years_between(a::Union{ZonedDateTime,DateTime,Date}, b::Union{ZonedDateTime,DateTime,Date})
    return (seconds_between(a, b) / (365.25 * 86400))
end


const global_base_date = ZonedDateTime(1970,1,1,0,0,0,tz"UTC")
"""
    period_length(a::Dates.DatePeriod, base::Date = global_base_date)

Period length is designed to convert `TimePeriod` objects to a float in a consistent way to `zdt2unix`.
So effectively the seconds_between method is calculated with start and end dates being those
at the start and end of a `Dates.DatePeriod`. This is slightly complicated because
a period like `Month(3)` might have slightly different numbers of total days depending
on when in the year it is. So a base date has to be input. The period is then
measured starting from this base date.

### Inputs
* `period` - A period.
* `base` - A date from which the period will be measured from.
### Returns
* A `Real`.
"""
function period_length(period::Dates.DatePeriod)
    return years_between(global_base_date+period, global_base_date)
end

"""
    years_from_global_base_date(a)

The number of years between a date and the global base date (1970-01-01 00:00:00 UTC).

This is used internally to convert dates to floats for use in `PE_Function` bases.

### Inputs
* `a` - A `ZonedDateTime`, `DateTime`, or `Date`

### Returns
* A `Real` representing the number of years since 1970-01-01.

### Example
```julia
years_from_global_base_date(Date(1970, 1, 1))  # 0.0
years_from_global_base_date(Date(2020, 1, 1))  # approximately 50.0
```
"""
years_from_global_base_date(a::Union{ZonedDateTime,DateTime,Date}) = years_between(a, global_base_date)
