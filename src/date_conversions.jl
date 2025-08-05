

function unix2zdt(seconds_since_1970::Integer, tzz = tz"UTC")
    if seconds_since_1970 > 3405888000
        error("This seems wrong unless you are far into the future. Are you sure you are not inputting microseconds since 1970?")
    end
    return ZonedDateTime(unix2datetime(seconds_since_1970), tzz)
end
function unix2dt(seconds_since_1970::Integer)
    if seconds_since_1970 > 3405888000
        error("This seems wrong unless you are far into the future. Are you sure you are not inputting microseconds since 1970?")
    end
    return unix2datetime(seconds_since_1970)
end
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

function days_between(a::Union{ZonedDateTime,DateTime,Date}, b::Union{ZonedDateTime,DateTime,Date})
    return (seconds_between(a, b) / 86400)
end

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

years_from_global_base_date(a::Union{ZonedDateTime,DateTime,Date}) = years_between(a, global_base_date)
