const tol = 10*eps()
const days_per_year = 365.2422
const global_base_date = Date(2000,1,1)
const global_base_date_as_day = convert(Dates.Day, global_base_date)

function years_between(a::Union{DateTime,Date}, b::Union{DateTime,Date})
    return (Dates.days(a) -Dates.days(b))/ days_per_year
end

function years_from_global_base(a::Union{DateTime,Date})
    return years_between(a, global_base_date)
end


"""
    Period length is designed to convert TimePeriod objects to a float in a consistent way to years_from_global_base
"""
function period_length(a::Dates.DatePeriod, base::Date = global_base_date)
    return years_between(base+a, base)
end
